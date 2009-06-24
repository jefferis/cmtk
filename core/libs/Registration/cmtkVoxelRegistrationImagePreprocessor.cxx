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

#include <cmtkVoxelRegistration.h>
#include <cmtkCommandLine.h>

namespace
cmtk
{

VoxelRegistration::ImagePreprocessor::ImagePreprocessor()
  : m_DataClassString( NULL ),
    m_PaddingFlag( false ),
    m_PaddingValue( 0 ),
    m_LowerThresholdActive( false ),
    m_LowerThresholdValue( 0 ),
    m_UpperThresholdActive( false ),
    m_UpperThresholdValue( 0 ),
    m_PruneHistogramBins( 0 ),
    m_CropIndex( NULL ),
    m_CropWorld( NULL ),
    m_AutoCropFlag( false ),
    m_AutoCropLevel( 0 )
{
}

void 
VoxelRegistration::ImagePreprocessor::AttachToCommandLine
( CommandLine& cl, const char* name, const char* key )
{
  char buffer[64];

  cl.BeginGroup( name, strcat( strcpy( buffer, name ), " Image Preprocessing" ) );
  
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "class-" ), key ) ), &this->m_DataClassString, "Data class: grey (default), binary, or label" );
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "pad-" ), key ) ), &this->m_PaddingValue, "Padding value", &this->m_PaddingFlag );
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "thresh-min-" ), key ) ), &this->m_LowerThresholdValue, "Minimum value truncation threshold", &this->m_LowerThresholdActive );
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "thresh-max-" ), key ) ), &this->m_UpperThresholdValue, "Maximum value truncation threshold", &this->m_UpperThresholdActive );
  
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "prune-histogram-" ), key ) ), &this->m_PruneHistogramBins, "Number of bins for histogram-based pruning [default: no pruning]" );
  
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "crop-index-" ), key ) ), &this->m_CropIndex, "Cropping region in pixel index coordinates [parsed as %d,%d,%d]" );
  cl.AddOption( CommandLine::Key( strcat( strcpy( buffer, "crop-world-" ), key ) ), &this->m_CropWorld, "Cropping region in world coordinates [parsed as %f,%f,%f]" );
  
  cl.EndGroup();
}

UniformVolume*
VoxelRegistration::ImagePreprocessor::GetProcessedImage( const UniformVolume* original )
{
  UniformVolume* volume = original->Clone();

  if ( this->m_DataClassString )
    {
    volume->GetData()->SetDataClass( StringToDataClass( this->m_DataClassString ) );
    }

  if ( this->m_PaddingFlag ) 
    {
    volume->GetData()->SetPaddingValue( this->m_PaddingValue );
    }

  if ( this->m_LowerThresholdActive || this->m_UpperThresholdActive )
    {
    volume->GetData()->Threshold( this->m_LowerThresholdValue, this->m_UpperThresholdValue );
    }

  if ( this->m_PruneHistogramBins )
    {
    volume->GetData()->PruneHistogram( true /*pruneHi*/, false /*pruneLo*/, this->m_PruneHistogramBins );
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

} // namespace cmtk
