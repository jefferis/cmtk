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

#include <cmtkImagePairSimilarityMeasure.h>

#include <cmtkUniformVolumeInterpolator.h>

#include <cmtkLinearInterpolator.h>
#include <cmtkNearestNeighborInterpolator.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairSimilarityMeasure::ImagePairSimilarityMeasure
( const UniformVolume::SmartPtr& refVolume, const UniformVolume::SmartPtr& fltVolume )
{
  this->m_ReferenceData = refVolume->GetData();
  this->m_FloatingData = fltVolume->GetData();

  switch ( this->m_FloatingData->GetDataClass() ) 
    {
    case DATACLASS_UNKNOWN :
    case DATACLASS_GREY :
      this->m_FloatingImageInterpolator = 
	cmtk::UniformVolumeInterpolatorBase::SmartPtr( new cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear>( fltVolume ) );
    case DATACLASS_LABEL :
      this->m_FloatingImageInterpolator = 
	cmtk::UniformVolumeInterpolatorBase::SmartPtr( new cmtk::UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor>( fltVolume ) );
    }
}

ImagePairSimilarityMeasure::ImagePairSimilarityMeasure
( const Self& other )
{
  this->m_ReferenceData = other.m_ReferenceData;
  this->m_FloatingData = other.m_FloatingData;
  this->m_FloatingImageInterpolator = other.m_FloatingImageInterpolator;
}

ImagePairSimilarityMeasure::ImagePairSimilarityMeasure
( Self& other, const bool copyData )
{
  if ( copyData )
    {
    }
  else
    {
    this->m_ReferenceData = other.m_ReferenceData;
    this->m_FloatingData = other.m_FloatingData;
    this->m_FloatingImageInterpolator = other.m_FloatingImageInterpolator;
    }
}

} // namespace cmtk
