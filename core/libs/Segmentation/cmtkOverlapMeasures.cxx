/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <Segmentation/cmtkOverlapMeasures.h>

#include <System/cmtkProgress.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

OverlapMeasures::OverlapMeasures( const std::vector<TypedArray::SmartPtr>& dataArrays )
{
  this->m_DataArrays = dataArrays;
  this->m_MaxLabelValue = 0;
  for ( size_t i = 0; i < this->m_DataArrays.size(); ++i )
    {
    const Types::DataItemRange range = this->m_DataArrays[i]->GetRange();
    this->m_MaxLabelValue = std::max<unsigned int>( this->m_MaxLabelValue, static_cast<unsigned int>( range.m_UpperBound ) );
    }
  
  // set size limits
  this->m_NumberOfImages = this->m_DataArrays.size();

  this->m_NumberOfPixels = this->m_DataArrays[0]->GetDataSize();
  for ( size_t i = 1; i < this->m_NumberOfImages; ++i )
    {
    this->m_NumberOfPixels = std::min( this->m_NumberOfPixels, this->m_DataArrays[i]->GetDataSize() );
    }
}

size_t
OverlapMeasures::ComputeGroupwiseOverlap
( const int firstLabel, const int numberOfLabels, double& overlapEqualWeighted, double& overlapVolumeWeighted, double& overlapInverseWeighted ) const
{
  // compute label volumes per image
  std::vector< std::vector< unsigned int > > volumeTable( numberOfLabels );
  for ( int label = 0; label < numberOfLabels; ++label )
    {
    volumeTable[label].resize( this->m_NumberOfImages, 0 );
    }

  std::vector<bool> labelExists( numberOfLabels );
  std::fill( labelExists.begin(), labelExists.end(), false );

  for ( size_t i = 0; i < this->m_NumberOfImages; ++i )
    {
    for ( size_t px = 0; px < this->m_NumberOfPixels; ++px )
      {
      Types::DataItem l;
      if ( this->m_DataArrays[i]->Get( l, px ) )
	{
	const int thisLabel = static_cast<int>(l) - firstLabel;
	if ( (thisLabel >= 0) && (thisLabel < numberOfLabels) )
	  {
	  ++volumeTable[thisLabel][i];
	  labelExists[thisLabel] = true;
	  }
	}
      }
    }

  size_t numberOfLabelsIncluded = 0;
  for ( int label = 0; label < numberOfLabels; ++label )
    {
    if ( labelExists[label] )
      ++numberOfLabelsIncluded;
    }
  if ( ! numberOfLabelsIncluded )
    return numberOfLabelsIncluded;

  // run overlap computation
  const size_t progressPixels = 100000;
  Progress::Begin( 0, this->m_NumberOfPixels, progressPixels, "Overlap computation" );

#ifdef _OPENMP
  const size_t numberOfThreads = omp_get_max_threads();
#else
  const size_t numberOfThreads = 1;
#endif

  std::vector<int> labelsPerPixel( numberOfThreads * this->m_NumberOfImages );
  
  const size_t sumsPerThread = numberOfLabels * this->m_NumberOfImages * (this->m_NumberOfImages-1) / 2;
  std::vector<double> sumsMin( numberOfThreads * sumsPerThread, 0.0 );
  std::vector<double> sumsMax( numberOfThreads * sumsPerThread, 0.0 );

#pragma omp parallel for  
  for ( int px = 0; px < static_cast<int>( this->m_NumberOfPixels ); ++px )
    {
    if ( (px % progressPixels) == 0 )
#ifdef _OPENMP
      if ( !omp_get_thread_num() )
#endif
	Progress::SetProgress( px );

#ifdef _OPENMP
    const size_t labelsPerPixelOfs = this->m_NumberOfImages * omp_get_thread_num();
#else
    const size_t labelsPerPixelOfs = 0;
#endif

    Types::DataItem l;
    for ( size_t i = 0; i < this->m_NumberOfImages; ++i )
      {
      labelsPerPixel[labelsPerPixelOfs+i] = -1;
      if ( this->m_DataArrays[i]->Get( l, px ) ) 
	{
	const int thisLabel = static_cast<int>(l) - firstLabel;
	if ( (thisLabel >= 0) && (thisLabel < numberOfLabels) )
	  labelsPerPixel[labelsPerPixelOfs+i] = thisLabel;
	}
      }
    
#ifdef _OPENMP
    size_t ofs = sumsPerThread * omp_get_thread_num();
#else
    size_t ofs = 0;
#endif
    for ( int label = 0; label < numberOfLabels; ++label )
      {
      if ( labelExists[label] )
	{
	for ( size_t i = 0; i < this->m_NumberOfImages; ++i )
	  {
	  const int wi = (labelsPerPixel[labelsPerPixelOfs+i] == label)?1:0;
	  for ( size_t j = 0; j < i; ++j, ++ofs )
	    {
	    const int wj = (labelsPerPixel[labelsPerPixelOfs+j] == label)?1:0;
	    
	    sumsMin[ofs] += std::min( wi, wj );
	    sumsMax[ofs] += std::max( wi, wj );
	    }
	  }
	}
      }
    }
  
#ifdef _OPENMP
  size_t dstOfs = sumsPerThread;
  for ( size_t thread = 1; thread < numberOfThreads; ++thread )
    {
    size_t thrOfs = 0;
    for ( size_t i = 0; i < sumsPerThread; ++i, ++thrOfs, ++dstOfs )
      {
      sumsMin[thrOfs] += sumsMin[dstOfs];
      sumsMax[thrOfs] += sumsMax[dstOfs];
      }
    }
#endif
  
  Progress::Done(); // rest should go fast
  
  double overlap_min_equal = 0, overlap_max_equal = 0;
  double overlap_min_volume = 0, overlap_max_volume = 0;
  double overlap_min_inverse = 0, overlap_max_inverse = 0;

  size_t ofs = 0;
  for ( int label = 0; label < numberOfLabels; ++label )
    {
    if ( labelExists[label] )
      {
      for ( size_t i = 0; i < this->m_NumberOfImages; ++i )
	{
	const unsigned int vi = volumeTable[label][i];
	for ( size_t j = 0; j < i; ++j, ++ofs )
	  {
	  overlap_min_volume += sumsMin[ofs];
	  overlap_max_volume += sumsMax[ofs];
	  
	  const unsigned int vij = vi + volumeTable[label][j];
	  if ( vij )
	    {
	    overlap_min_equal += sumsMin[ofs] / vij;
	    overlap_max_equal += sumsMax[ofs] / vij;
	    
	    overlap_min_inverse += (sumsMin[ofs] / vij) / vij;
	    overlap_max_inverse += (sumsMax[ofs] / vij) / vij;
	    }
	  }
	}
      }
    }
  
  if ( overlap_max_equal )
    overlapEqualWeighted = overlap_min_equal / overlap_max_equal;
  if ( overlap_max_volume )
    overlapVolumeWeighted = overlap_min_volume / overlap_max_volume;
  if ( overlap_max_inverse )
    overlapInverseWeighted = overlap_min_inverse / overlap_max_inverse;

  return numberOfLabelsIncluded;
}

double
OverlapMeasures::ComputePairwiseOverlapMinMax
( double& overlap_min, double& overlap_max, const TypedArray::SmartPtr& data0, const TypedArray::SmartPtr& data1, const int label ) const
{
  overlap_min = overlap_max = 0;
  for ( size_t i = 0; i < this->m_NumberOfPixels; ++i )
    {
    Types::DataItem l0, l1;
    if ( ! data0->Get( l0, i ) ) l0 = -1;
    if ( ! data1->Get( l1, i ) ) l1 = -1;

    const int w0 = (l0 == label)?1:0;
    const int w1 = (l1 == label)?1:0;

    overlap_min += std::min( w0, w1 );
    overlap_max += std::max( w0, w1 );
    }
  return 0;
}

} // namespace cmtk

