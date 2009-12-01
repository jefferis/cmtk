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

#include <cmtkImagePairSimilarityMeasureMSD.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairSimilarityMeasureMSD::ImagePairSimilarityMeasureMSD
( const UniformVolume::SmartPtr& refVolume, const UniformVolume::SmartPtr& fltVolume ) :
  ImagePairSimilarityMeasure( refVolume, fltVolume )
{}

ImagePairSimilarityMeasureMSD::ImagePairSimilarityMeasureMSD
( Self& other, const bool copyData ) :
  ImagePairSimilarityMeasure( other, copyData )
{
  this->m_SumOfDifferences = other.m_SumOfDifferences;
  this->m_NumberOfSamples = other.m_NumberOfSamples;
}

ImagePairSimilarityMeasureMSD::ImagePairSimilarityMeasureMSD
( const Self& other ) :
  ImagePairSimilarityMeasure( other )
{
  this->m_SumOfDifferences = other.m_SumOfDifferences;
  this->m_NumberOfSamples = other.m_NumberOfSamples;
}

void
ImagePairSimilarityMeasureMSD::CopyUnsafe
( const ImagePairSimilarityMeasureMSD& other, const bool )
{
  this->m_SumOfDifferences = other.m_SumOfDifferences;
  this->m_NumberOfSamples = other.m_NumberOfSamples;
}

} // namespace cmtk
