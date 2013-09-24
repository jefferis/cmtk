/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
//
//  This file is part of the Computational Morphometry Toolkit.
//
//  http://www.nitrc.org/projects/cmtk/
//
//  The Computational Morphometry Toolkit is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  The Computational Morphometry Toolkit is distributed in the hope that it
//  will be useful, but WITHOUT ANY WARRANTY; without even the implied
//  warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with the Computational Morphometry Toolkit.  If not, see
//  <http://www.gnu.org/licenses/>.
//
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <Registration/cmtkVoxelRegistration.h>

#include <System/cmtkCommandLine.h>

#include <Base/cmtkTypedArrayFunctionHistogramEqualization.h>
#include <Base/cmtkDataGridFilter.h>

namespace
cmtk
{

VoxelRegistration::ImagePreprocessor::ImagePreprocessor( const std::string& name, const std::string& key )
  : m_DataClassString( NULL ),
    m_DataClass( DATACLASS_GREY ),
    m_PaddingFlag( false ),
    m_PaddingValue( 0 ),
    m_LowerThresholdActive( false ),
    m_LowerThresholdValue( -CMTK_ITEM_MAX ),
    m_UpperThresholdActive( false ),
    m_UpperThresholdValue( CMTK_ITEM_MAX ),
    m_UsePruneHistogramBins( false ),
    m_PruneHistogramBins( 0 ),
    m_HistogramEqualization( false ),
    m_SobelFilter( false ),
    m_CropIndex( NULL ),
    m_CropWorld( NULL ),
    m_AutoCropFlag( false ),
    m_AutoCropLevel( 0 ),
    m_Name( name ),
    m_Key( key )
{
}

void 
VoxelRegistration::ImagePreprocessor::AttachToCommandLine
( CommandLine& cl )
{
  cl.BeginGroup( this->m_Name, this->m_Name + " Image Preprocessing" )->SetProperties( CommandLine::PROPS_NOXML );
  
  cl.AddOption( CommandLine::Key( std::string( "class-" ) + this->m_Key ), &this->m_DataClassString, "Data class: grey (default) or label" );
  cl.AddOption( CommandLine::Key( std::string( "pad-" ) + this->m_Key ), &this->m_PaddingValue, "Padding value", &this->m_PaddingFlag );

  cl.AddOption( CommandLine::Key( std::string( "thresh-min-" ) + this->m_Key ), &this->m_LowerThresholdValue, "Minimum value truncation threshold", &this->m_LowerThresholdActive );
  cl.AddOption( CommandLine::Key( std::string( "thresh-max-" ) + this->m_Key ), &this->m_UpperThresholdValue, "Maximum value truncation threshold", &this->m_UpperThresholdActive );
  
  cl.AddOption( CommandLine::Key( std::string( "prune-histogram-" ) + this->m_Key ), &this->m_PruneHistogramBins, "Number of bins for histogram-based pruning", &this->m_UsePruneHistogramBins );
  cl.AddSwitch( CommandLine::Key( std::string( "histogram-equalization-" ) + this->m_Key ), &this->m_HistogramEqualization, true, "Apply histogram equalization" );
  cl.AddSwitch( CommandLine::Key( std::string( "sobel-filter-" ) + this->m_Key ), &this->m_SobelFilter, true, "Apply Sobel edge detection filter" );
  
  cl.AddOption( CommandLine::Key( std::string( "crop-index-" ) + this->m_Key ), &this->m_CropIndex, "Cropping region in pixel index coordinates [parsed as %d,%d,%d,%d,%d,%d for i0,j0,k0,i1,j1,k1]" );
  cl.AddOption( CommandLine::Key( std::string( "crop-world-" ) + this->m_Key ), &this->m_CropWorld, "Cropping region in world coordinates [parsed as %f,%f,%f,%f,%f,%f for x0,y0,z0,x1,y1,z1]" );
  
  cl.EndGroup();
}

UniformVolume::SmartPtr
VoxelRegistration::ImagePreprocessor::GetProcessedImage( const UniformVolume* original )
{
  UniformVolume::SmartPtr volume( original->Clone() );
  TypedArray::SmartPtr data = volume->GetData();

  if ( this->m_DataClassString )
    {
    this->m_DataClass = StringToDataClass( this->m_DataClassString );
    data->SetDataClass( this->m_DataClass );
    }

  if ( this->m_PaddingFlag ) 
    {
    data->SetPaddingValue( this->m_PaddingValue );
    }

  if ( this->m_LowerThresholdActive || this->m_UpperThresholdActive )
    {
    data->Threshold( Types::DataItemRange( this->m_LowerThresholdValue, this->m_UpperThresholdValue ) );
    }

  if ( this->m_PruneHistogramBins )
    {
    data->PruneHistogram( true /*pruneHi*/, false /*pruneLo*/, this->m_PruneHistogramBins );
    }
  
  if ( this->m_HistogramEqualization ) 
    {
    data->ApplyFunctionObject( TypedArrayFunctionHistogramEqualization( *data ) );
    }

  if ( this->m_SobelFilter ) 
    {
    volume->SetData( DataGridFilter( volume ).GetDataSobelFiltered() );
    }

  if ( this->m_CropIndex )
    {
    int cropFrom[3], cropTo[3];
    if ( 6 != sscanf( this->m_CropIndex, "%6d,%6d,%6d,%6d,%6d,%6d", cropFrom, cropFrom+1, cropFrom+2, cropTo, cropTo+1, cropTo+2 ) ) 
      {
      StdErr << "Option index coordinate cropping expects six integer parameters but got '" << this->m_CropIndex << "'\n";
      exit( 1 );
      }
    
    for ( int dim=0; dim<3; ++dim ) 
      {
      if ( cropTo[dim] < 0 ) 
	{
	cropTo[dim] = volume->GetDims()[dim] + cropTo[dim] + 1;
	}
      }
    volume->CropRegion() = DataGrid::RegionType( DataGrid::IndexType::FromPointer( cropFrom ), DataGrid::IndexType::FromPointer( cropTo ) );
    }
  
  if ( this->m_CropWorld )
    {
    float crop[6];
    if ( 6 != sscanf( this->m_CropWorld, "%15f,%15f,%15f,%15f,%15f,%15f", crop, crop+1, crop+2, crop+3, crop+4, crop+5 ) ) 
      {
      StdErr << "Option world coordinate cropping expects six floating-point parameters but got '" << this->m_CropWorld << "'\n";
      exit( 1 );
      }
    
    Types::Coordinate realCropFrom[3], realCropTo[3];
    for ( int dim=0; dim<3; ++dim ) 
      {
      realCropFrom[dim] = crop[dim];
      if ( crop[3+dim] < 0 ) 
	{
	realCropTo[dim] = volume->m_Size[dim] + crop[3+dim];
	}
      else
	{
	realCropTo[dim] = crop[3+dim];
	}
      }
    volume->SetHighResCropRegion( UniformVolume::CoordinateRegionType( UniformVolume::CoordinateRegionType::IndexType::FromPointer( realCropFrom ), UniformVolume::CoordinateRegionType::IndexType::FromPointer( realCropTo ) ) );
    }

  if ( this->m_AutoCropFlag )
    {
    volume->AutoCrop( this->m_AutoCropLevel, true /*recrop*/ );
    }

  return volume;
}

void 
VoxelRegistration::ImagePreprocessor::WriteSettings
( ClassStreamOutput& stream ) const
{
  stream.Begin( std::string( "preprocessing_" ) + this->m_Key );
  switch ( this->m_DataClass ) 
    {
    case DATACLASS_GREY:
      stream.WriteString( "dataclass", "GreyLevel" );
      break;
    case DATACLASS_LABEL:
      stream.WriteString( "dataclass", "LabelField" );
      break;
    default:
      stream.WriteString( "dataclass", "Unknown" );
      break;
    }

  if ( this->m_PaddingFlag )
    stream.WriteDouble( "padding_value", this->m_PaddingValue );

  if ( this->m_LowerThresholdActive )
    stream.WriteDouble( "thresh_lower", this->m_LowerThresholdValue );    
    
  if ( this->m_UpperThresholdActive )
    stream.WriteDouble( "thresh_upper", this->m_UpperThresholdValue );    
    
  if ( this->m_PruneHistogramBins )
    stream.WriteInt( "prune_histogram_bins", this->m_PruneHistogramBins );

  if ( this->m_HistogramEqualization ) 
    stream.WriteBool( "histogram_equalization", true );

  if ( this->m_SobelFilter ) 
    stream.WriteBool( "sobel_filter", true );

  if ( this->m_CropIndex )
    stream.WriteString( "crop_index", this->m_CropIndex );

  if ( this->m_CropWorld )
    stream.WriteString( "crop_world", this->m_CropWorld );

  if ( this->m_AutoCropFlag )
    stream.WriteDouble( "auto_crop_level", this->m_AutoCropLevel );

  stream.End();
}  

} // namespace cmtk
