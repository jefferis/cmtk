/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkImagePairRegistration.h>

#include <cmtkCommandLine.h>
#include <cmtkTypedArrayFunctionHistogramEqualization.h>

namespace
cmtk
{

ImagePairRegistration::ImagePreprocessor::ImagePreprocessor( const char* name, const char* key )
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
ImagePairRegistration::ImagePreprocessor::AttachToCommandLine
( CommandLine& cl )
{
  char buffer[64];
  cl.BeginGroup( this->m_Name, strcat( strcpy( buffer, this->m_Name ), " Image Preprocessing" ) )->SetProperties( CommandLine::PROPS_NOXML );
  
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "class-" ), this->m_Key ) ), &this->m_DataClassString, "Data class: grey (default) or label" );
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "pad-" ), this->m_Key ) ), &this->m_PaddingValue, "Padding value", &this->m_PaddingFlag );

  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "thresh-min-" ), this->m_Key ) ), &this->m_LowerThresholdValue, "Minimum value truncation threshold", &this->m_LowerThresholdActive );
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "thresh-max-" ), this->m_Key ) ), &this->m_UpperThresholdValue, "Maximum value truncation threshold", &this->m_UpperThresholdActive );
  
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "prune-histogram-" ), this->m_Key ) ), &this->m_PruneHistogramBins, "Number of bins for histogram-based pruning", &this->m_UsePruneHistogramBins );
  cl.AddSwitch( CommandLine::Key( strcat( strcpy( buffer, "histogram-equalization-" ), this->m_Key ) ), &this->m_HistogramEqualization, true, "Apply histogram equalization" );
  cl.AddSwitch( CommandLine::Key( strcat( strcpy( buffer, "sobel-filter-" ), this->m_Key ) ), &this->m_SobelFilter, true, "Apply Sobel edge detection filter" );
  
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "crop-index-" ), this->m_Key ) ), &this->m_CropIndex, "Cropping region in pixel index coordinates [parsed as %d,%d,%d,%d,%d,%d for i0,j0,k0,i1,j1,k1]" );
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "crop-world-" ), this->m_Key ) ), &this->m_CropWorld, "Cropping region in world coordinates [parsed as %f,%f,%f,%f,%f,%f for x0,y0,z0,x1,y1,z1]" );
  
  cl.EndGroup();
}

UniformVolume*
ImagePairRegistration::ImagePreprocessor::GetProcessedImage( const UniformVolume* original )
{
  UniformVolume* volume = original->Clone();
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
    data->Threshold( this->m_LowerThresholdValue, this->m_UpperThresholdValue );
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
    volume->ApplySobelFilter();
    }

  if ( this->m_CropIndex )
    {
    int crop[6];
    if ( 6 != sscanf( this->m_CropIndex, "%d,%d,%d,%d,%d,%d", crop, crop+1, crop+2, crop+3, crop+4, crop+5 ) ) 
      {
      StdErr << "Option index coordinate cropping expects six integer parameters but got '" << this->m_CropIndex << "'\n";
      exit( 1 );
      }

    for ( int dim=0; dim<3; ++dim ) 
      {
      if ( crop[3+dim] < 0 ) 
	{
	crop[3+dim] = volume->GetDims()[dim] + crop[3+dim] + 1;
	}
      }
    volume->SetCropRegion( crop, crop+3 );
    }
  
  if ( this->m_CropWorld )
    {
    float crop[6];
    if ( 6 != sscanf( this->m_CropWorld, "%f,%f,%f,%f,%f,%f", crop, crop+1, crop+2, crop+3, crop+4, crop+5 ) ) 
      {
      StdErr << "Option world coordinate cropping expects six floating-point parameters but got '" << this->m_CropWorld << "'\n";
      exit( 1 );
      }
    
    Types::Coordinate realCrop[6];
    for ( int dim=0; dim<3; ++dim ) 
      {
      realCrop[dim] = crop[dim];
      if ( crop[3+dim] < 0 ) 
	{
	realCrop[3+dim] = volume->Size[dim] + crop[3+dim];
	}
      else
	{
	realCrop[3+dim] = crop[3+dim];
	}
      }
    volume->SetCropRegion( realCrop, realCrop+3 );
    }

  if ( this->m_AutoCropFlag )
    {
    volume->AutoCrop( this->m_AutoCropLevel, true /*recrop*/ );
    }

  return volume;
}

void 
ImagePairRegistration::ImagePreprocessor::WriteSettings
( ClassStream& stream ) const
{
  char buffer[64];
  stream.Begin( strcat( strcpy( buffer, "preprocessing_" ), this->m_Key ) );
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
