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

namespace
cmtk
{

VoxelRegistration::ImagePreprocessor::ImagePreprocessor()
  : m_DataClass( DATACLASS_GREY ),
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

UniformVolume*
VoxelRegistration::ImagePreprocessor::GetProcessedImage( const UniformVolume* original )
{
  UniformVolume* volume = original->Clone();

  volume->GetData()->SetDataClass( this->m_DataClass );
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
