/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include <cmtkReformatVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairSimilarityMeasure::ImagePairSimilarityMeasure
( const UniformVolume::SmartConstPtr& refVolume, const UniformVolume::SmartConstPtr& fltVolume, const Interpolators::InterpolationEnum interpolation )
  : m_InterpolationMethod( interpolation )
{
  this->SetReferenceVolume( refVolume );
  this->SetFloatingVolume( fltVolume );
}

void
ImagePairSimilarityMeasure::SetReferenceVolume( const UniformVolume::SmartConstPtr& refVolume )
{
  this->m_ReferenceVolume = refVolume;
  this->m_ReferenceData = this->m_ReferenceVolume->GetData();
}

void
ImagePairSimilarityMeasure::SetFloatingVolume( const UniformVolume::SmartConstPtr& fltVolume )
{
  this->m_FloatingVolume = fltVolume;
  this->m_FloatingData = fltVolume->GetData();
  
  if ( this->m_InterpolationMethod == Interpolators::DEFAULT )
    {
    // decide based on floating image data class.
    switch ( this->m_FloatingData->GetDataClass() ) 
      {
      case DATACLASS_UNKNOWN :
      case DATACLASS_GREY :
	this->m_InterpolationMethod = Interpolators::LINEAR;
	this->m_FloatingImageInterpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear>( fltVolume ) );
	break;
      case DATACLASS_LABEL :
	this->m_InterpolationMethod = Interpolators::NEAREST_NEIGHBOR;
	this->m_FloatingImageInterpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new cmtk::UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor>( fltVolume ) );
	break;
      }
    }
  else
    {
    this->m_FloatingImageInterpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( ReformatVolume::CreateInterpolator( this->m_InterpolationMethod, fltVolume ) );
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
    StdErr << "Not implemented: " << __FILE__ << ":" << __LINE__ << "\n";
    exit(1);
    }
  else
    {
    this->m_ReferenceData = other.m_ReferenceData;
    this->m_FloatingData = other.m_FloatingData;
    this->m_FloatingImageInterpolator = other.m_FloatingImageInterpolator;
    }
}

} // namespace cmtk
