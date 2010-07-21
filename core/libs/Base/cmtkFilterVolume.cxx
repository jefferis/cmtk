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

#include "cmtkFilterVolume.h"

#include "Base/cmtkMathUtil.h"
#include "Base/cmtkFilterMask.h"
#include "Base/cmtkHistogram.h"
#include "Base/cmtkJointHistogram.h"

#include "System/cmtkProgress.h"
#include "System/cmtkThreads.h"
#include "System/cmtkException.h"

#ifdef _OPENMP
#  include <omp.h>
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{

TypedArray::SmartPtr
FilterVolume::GaussianFilter
( const UniformVolume* volume, const Units::GaussianSigma& kernelWidth, const Types::Coordinate radius, const TypedArray* maskData )
{
  const TypedArray* inputData = volume->GetData();
  if ( ! inputData ) 
    throw( Exception( "Missing image data" ) );
  
  TypedArray::SmartPtr filtered = TypedArray::Create( inputData->GetType(), inputData->GetDataSize() );
  
  const DataGrid::IndexType& dims = volume->m_Dims;
  FilterMask<3> filter( dims, volume->GetDelta(), radius, FilterMask<3>::Gaussian( kernelWidth ) );

  const int dimsX = dims[AXIS_X];
  const int dimsY = dims[AXIS_Y];
  const int dimsZ = dims[AXIS_Z];
  
  Progress::Begin( 0, dimsZ, 1, "Gaussian Filter" );

#pragma omp parallel for
  for ( int z = 0; z < dimsZ; ++z )
    {
    size_t offset = z * dimsX * dimsY;
    Progress::SetProgress( z );
    for ( int y = 0; y < dimsY; ++y )
      for ( int x = 0; x < dimsX; ++x, ++offset ) 
	{
	Types::DataItem average = 0.0, weight = 0.0;
	  
	Types::DataItem maskValue = 0.0;
	if ( maskData )
	  {
	  maskData->Get( maskValue, offset );
	  }
	else
	  {
	  maskValue = 1.0;
	  }
	  
	if ( maskValue ) 
	  {
	  FilterMask<3>::const_iterator it = filter.begin();
	  while ( it != filter.end() ) 
	    {
	    const int xx = x + it->Location[0];
	    const int yy = y + it->Location[1];
	    const int zz = z + it->Location[2];
		    
	    if ( (xx>=0) && (yy>=0) && (zz>=0) && (xx < (int)dimsX) && (yy < (int)dimsY) && (zz < (int)dimsZ) ) 
	      {
	      const size_t srcOffset = volume->GetOffsetFromIndex( xx, yy, zz );
	      Types::DataItem value;
	      if ( inputData->Get( value, srcOffset ) ) 
		{
		average += it->Coefficient * value;
		weight += it->Coefficient;
		}
	      }
	    ++it;	    
	    }
	  }
	if ( weight > 0.0 ) 
	  {
	  filtered->Set( average / weight, offset );
	  } 
	else 
	  {
	  filtered->SetPaddingAt( offset );
	  }
	}
    }
  
  Progress::Done();

  return filtered;
}

void 
printBlock
( Types::DataItem block[COUPE_BLOCK_SIZE] )
{
  for ( int i = 0; i < COUPE_BLOCK_SIZE; i++ )
    std::cout << block[i] << "   ";
}

TypedArray::SmartPtr
FilterVolume
::RohlfingFilter
( const UniformVolume* volume, const TypedArray* subjectData,
  const TypedArray* maskData, const Units::GaussianSigma& iFilterSigma,
  const Units::GaussianSigma& filterWidth, const Types::Coordinate filterRadius )
{
  const TypedArray* inputData = volume->GetData();
  if ( ! inputData )
    throw( Exception( "Missing image data" ) );
 
  const Types::DataItemRange rangeSubj = subjectData->GetRange();
  
  const size_t numBins = 1024;
#ifdef _OPENMP
  const size_t maxThreads = omp_get_max_threads();
  std::vector<Histogram<Types::DataItem>::SmartPtr> histograms( maxThreads );
  for ( size_t thread = 0; thread < maxThreads; ++thread )
    {
    histograms[thread] = Histogram<Types::DataItem>::SmartPtr( new Histogram<Types::DataItem>( numBins ) );
    histograms[thread]->SetRange( rangeSubj );
    }
#else // #ifdef _OPENMP
  Histogram<Types::DataItem> histogram( numBins );
  histogram.SetRange( rangeSubj );
#endif // #ifdef _OPENMP

  const size_t iKernelRadius = 1 + static_cast<size_t>( 2 * iFilterSigma.Value() * numBins );
  Types::DataItem* iKernel = Memory::AllocateArray<Types::DataItem>( iKernelRadius );
  if ( iKernelRadius > 1 )
    {
    const Types::DataItem normFactor = static_cast<Types::DataItem>( 1.0/(sqrt(2*M_PI) * iFilterSigma.Value() * numBins) ); // not really necessary since we normalize during convolution
    for ( size_t i = 0; i < iKernelRadius; ++i )
      {
      iKernel[i] = static_cast<Types::DataItem>( normFactor * exp( -MathUtil::Square( 1.0 * i / (iFilterSigma.Value()*numBins) ) / 2 ) );
      }
    }
  else
    {
    iKernel[0] = 1.0;
    }
  
  TypedArray::SmartPtr filtered = TypedArray::Create( inputData->GetType(), inputData->GetDataSize() );
  
  const DataGrid::IndexType& dims = volume->GetDims();
  FilterMask<3> filter( dims, volume->GetDelta(), filterRadius, FilterMask<3>::Gaussian( filterWidth ) );
  
  const unsigned int dimsX = dims[AXIS_X];
  const unsigned int dimsY = dims[AXIS_Y];
  const unsigned int dimsZ = dims[AXIS_Z];
  
  Progress::Begin( 0, dimsZ, 1, "Rohlfing Intensity-Consistent Filter" );

#pragma omp parallel for
  for ( unsigned int z = 0; z < dimsZ; ++z ) 
    {      
    size_t offset = z * dimsX * dimsY;

#ifdef _OPENMP
    const size_t threadIdx = omp_get_thread_num();
    Histogram<Types::DataItem>& histogram = *(histograms[threadIdx]);
    if ( ! threadIdx )
#endif // #ifdef _OPENMP
      Progress::SetProgress( z );

    for ( unsigned int y = 0; y < dimsY; ++y )
      for ( unsigned int x = 0; x < dimsX; ++x, ++offset ) 
	{
	Types::DataItem average = 0.0, weight = 0.0;
	    
	Types::DataItem maskValue = 1.0;
	if ( maskData )
	  maskData->Get( maskValue, offset );
	    
	Types::DataItem valueSubj;
	if ( maskValue && subjectData->Get( valueSubj, offset ) ) 
	  {
	  histogram.Reset();
	  histogram.AddWeightedSymmetricKernel( histogram.ValueToBin( valueSubj ), iKernelRadius, iKernel );
	  
	  for (  FilterMask<3>::const_iterator it = filter.begin(); it != filter.end(); ++it ) 
	    {
	    const int xx = x + it->Location[0];
	    const int yy = y + it->Location[1];
	    const int zz = z + it->Location[2];
		    
	    if ( (xx>=0) && (yy>=0) && (zz>=0) && (xx < (int)dimsX) && (yy < (int)dimsY) && (zz < (int)dimsZ) ) 
	      {
	      Types::DataItem value;
	      const size_t srcOffset = it->RelativeIndex + offset;
	      if ( inputData->Get( value, srcOffset ) ) 
		{			    
		Types::DataItem valueSubj;
		if ( subjectData->Get( valueSubj, srcOffset ) )
		  {
		  const size_t bin = histogram.ValueToBin( valueSubj );
		  const Types::DataItem prob = it->Coefficient * histogram[bin];
		  
		  average += value * prob;
		  weight += prob;
		  }
		}
	      }
	    }
	  }
	
	if ( weight > 0.0 )
	  {
	  filtered->Set( average / weight, offset );
	  } 
	else
	  {
	  filtered->SetPaddingAt( offset );
	  }
	}
    }
  
  Progress::Done();

  delete[] iKernel;
  
  return filtered;
}

TypedArray::SmartPtr
FilterVolume::StudholmeFilter
( const UniformVolume* volume, const TypedArray* subjectData,
  const TypedArray* averageData, const TypedArray* maskData,
  std::list<TypedArray::SmartPtr> imgList, const Types::DataItem binWidth, 
  const Units::GaussianSigma& filterWidth, const Types::Coordinate filterRadius )
{
  const TypedArray* inputData = volume->GetData();
  if ( ! inputData )
    throw( Exception( "Missing image data" ) );
 
  const Types::DataItemRange range = averageData->GetRange();
  const size_t numBins = std::min( 128, 1 + static_cast<int>((range.Width()) / binWidth) );

  TypedArray::SmartPtr filtered = TypedArray::Create( inputData->GetType(), inputData->GetDataSize() );
  
  const DataGrid::IndexType& dims = volume->GetDims();
  const unsigned int dimsX = dims[AXIS_X];
  const unsigned int dimsY = dims[AXIS_Y];
  const unsigned int dimsZ = dims[AXIS_Z];
  const unsigned int numberOfRows = dimsY * dimsZ;
  
  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  std::vector< JointHistogram<Types::DataItem> > histogramByThread( numberOfThreads );
  std::vector<FilterMask<3>::SmartPtr> filterByThread( numberOfThreads );
  
  for ( size_t idx = 0; idx < numberOfThreads; ++idx )
    {
    histogramByThread[idx].SetNumBins( numBins, numBins );
    histogramByThread[idx].SetRangeX( range );
    histogramByThread[idx].SetRangeY( range );
    
    FilterMask<3>::SmartPtr filter( new FilterMask<3>( dims, volume->GetDelta(), filterRadius, FilterMask<3>::Gaussian( filterWidth ) ) );    
    filterByThread[idx] = filter;
    }
  
  Progress::Begin( 0, numberOfRows, 1, "Studholme Intensity-Consistent Filter" );
#pragma omp parallel for
  for ( unsigned int row = 0; row < numberOfRows; ++row ) 
    {
    const unsigned int y = row % dimsY;
    const unsigned int z = row / dimsY;
    
    Progress::SetProgress( z );
    size_t offset = row * dimsX;
     
#ifdef _OPENMP
    const int thread = omp_get_thread_num();
#else
    const int thread = 0;
#endif

    JointHistogram<Types::DataItem>& histogram = histogramByThread[thread];
    FilterMask<3>& filter = *(filterByThread[thread]);
 
    for ( unsigned int x = 0; x < dimsX; ++x, ++offset ) 
      {
      Types::DataItem average = 0.0, weight = 0.0;
      histogram.Reset();
	    
      Types::DataItem maskValue = 1.0;
      if ( maskData )
	maskData->Get( maskValue, offset );
	    
      Types::DataItem valueAvg;
      if ( maskValue && averageData->Get( valueAvg, offset ) ) 
	{
		
	// first iteration over filter: compute consistency histogram
	FilterMask<3>::iterator it = filter.begin();
	for ( ; it != filter.end(); ++it ) 
	  {
	  const int xx = x + it->Location[0];
	  const int yy = y + it->Location[1];
	  const int zz = z + it->Location[2];
		    
	  if ( (xx>=0) && (yy>=0) && (zz>=0) && (xx < (int)dimsX) && (yy < (int)dimsY) && (zz < (int)dimsZ) ) 
	    {
	    it->Valid = true;	      

	    const size_t srcOffset = it->RelativeIndex + offset;
	    it->PixelIndex = srcOffset;
			
	    Types::DataItem valueAvgSrc, valueSubj;
	    if ( averageData->Get( valueAvgSrc, srcOffset ) ) 
	      {			    
	      const size_t binAvg = histogram.ValueToBinX( valueAvgSrc );
	      for ( std::list<TypedArray::SmartPtr>::iterator itImg = imgList.begin(); itImg != imgList.end(); ++itImg ) 
		{
		if ( (*itImg)->Get( valueSubj, srcOffset ) ) 
		  {
		  histogram.Increment( binAvg, histogram.ValueToBinY( valueSubj ) );
		  }
		}
		
	      }
	    }
	  else
	    {
	    it->Valid = false;
	    }		    
	  }
	  
	const size_t binX = histogram.ValueToBinX( valueAvg );
	const Types::DataItem avgHistogramValueInv = static_cast<Types::DataItem>( 1.0/histogram.ProjectToX( binX ) );
	  
	for ( it = filter.begin(); it != filter.end(); ++it ) 
	  {
	  if ( it->Valid ) 
	    {
	    Types::DataItem value;
	    if ( inputData->Get( value, it->PixelIndex ) ) 
	      {			    
	      Types::DataItem valueSubj;
	      if ( subjectData->Get( valueSubj, it->PixelIndex ) )
		{
		const size_t binY = histogram.ValueToBinY( valueSubj );
		  
		const Types::DataItem prob = static_cast<Types::DataItem>( it->Coefficient * avgHistogramValueInv * histogram.GetBin( binX, binY ) );
		  
		average += value * prob;
		weight += prob;
		}
	      }
	    }
	  }
	}
	
      if ( weight > 0.0 )
	{
	filtered->Set( average / weight, offset );
	} 
      else
	{
	filtered->SetPaddingAt( offset );
	}
      }
    }
  
  Progress::Done();
  
  return filtered;
}

TypedArray::SmartPtr
FilterVolume::StudholmeFilter
( const UniformVolume* volume, 
  std::list<TypedArray::SmartPtr> subjectData,
  const TypedArray* averageData, const TypedArray* maskData,
  std::list<TypedArray::SmartPtr> imgList, const Types::DataItem binWidth, 
  const Units::GaussianSigma& filterWidth, const Types::Coordinate filterRadius )
{
  const TypedArray* inputData = volume->GetData();
  if ( ! inputData ) 
    throw( Exception( "Missing image data" ) );
 
  const Types::DataItemRange range = averageData->GetRange();

  const size_t numBins = std::min( 128, 1 + static_cast<int>( range.Width() / binWidth ) );
  JointHistogram<Types::DataItem> histogram( numBins, numBins );
  histogram.SetRangeX( range );
  histogram.SetRangeY( range );
 
  TypedArray::SmartPtr filtered = TypedArray::Create( inputData->GetType(), inputData->GetDataSize() );
  
  const DataGrid::IndexType& dims = volume->GetDims();
  FilterMask<3> filter( dims, volume->GetDelta(), filterRadius, FilterMask<3>::Gaussian( filterWidth ) );
  
  const unsigned int dimsX = dims[AXIS_X];
  const unsigned int dimsY = dims[AXIS_Y];
  const unsigned int dimsZ = dims[AXIS_Z];

  Progress::Begin( 0, dimsZ, 1, "Studholme Intensity-Consistent Filter" );

  size_t offset = 0;
  for ( unsigned int z = 0; z < dimsZ; ++z ) 
    {
    Progress::SetProgress( z );
      
    for ( unsigned int y = 0; y < dimsY; ++y )
      for ( unsigned int x = 0; x < dimsX; ++x, ++offset ) 
	{
	Types::DataItem average = 0.0, weight = 0.0;
	histogram.Reset();
	    
	Types::DataItem maskValue = 1.0;
	if ( maskData )
	  maskData->Get( maskValue, offset );
	    
	Types::DataItem valueAvg;
	if ( maskValue && averageData->Get( valueAvg, offset ) ) 
	  {
	  // first iteration over filter: compute consistency histogram
	  for ( FilterMask<3>::iterator it = filter.begin(); it != filter.end(); ++it ) 
	    {
	    const unsigned int xx = x + it->Location[0];
	    const unsigned int yy = y + it->Location[1];
	    const unsigned int zz = z + it->Location[2];
		    
	    if ( (xx < dimsX) && (yy < dimsY) && (zz < dimsZ) ) 
	      {
	      it->Valid = true;	      
	      // since xx, yy, zz are unsigned, we need not check 
	      // for >= 0; this is taken care of by overflow (we 
	      // hope ;)
	      const size_t srcOffset = volume->GetOffsetFromIndex( xx, yy, zz );
	      Types::DataItem valueAvgSrc, valueSubj;
	      if ( averageData->Get( valueAvgSrc, srcOffset ) ) 
		{
		const size_t binAvg = histogram.ValueToBinX( valueAvgSrc );
		for ( std::list<TypedArray::SmartPtr>::iterator itImg = imgList.begin(); itImg != imgList.end(); ++itImg ) 
		  {
		  if ( (*itImg)->Get( valueSubj, srcOffset ) ) 
		    {
		    histogram.Increment( binAvg, histogram.ValueToBinY( valueSubj ) );
		    }
		  }		
		}
	      }
	    }
	  
	  Histogram<Types::DataItem>* avgHistogram = histogram.GetMarginalX();
	  const size_t binX = histogram.ValueToBinX( valueAvg );
	  
	  for ( FilterMask<3>::iterator it = filter.begin(); it != filter.end(); ++it ) 
	    {
	    const unsigned int xx = x + it->Location[0];
	    const unsigned int yy = y + it->Location[1];
	    const unsigned int zz = z + it->Location[2];
	    
	    if ( it->Valid ) 
	      {
	      it->Valid = false;
	      // since xx, yy, zz are unsigned, we need not check for
	      // >= 0; this is taken care of by overflow (we hope ;)
	      const size_t srcOffset = volume->GetOffsetFromIndex( xx, yy, zz );
	      Types::DataItem value;
	      if ( inputData->Get( value, srcOffset ) ) 
		{
		float prob = it->Coefficient;
		
		std::list<TypedArray::SmartPtr>::iterator subjectIt = subjectData.begin();
		while ( subjectIt != subjectData.end() ) 
		  {
		  Types::DataItem valueSubj;
		  if ( (*subjectIt)->Get( valueSubj, srcOffset ) )
		    {
		    const size_t binY = histogram.ValueToBinY( valueSubj );
		    prob *= histogram.GetBin( binX, binY ) / (*avgHistogram)[binX];
		    }
		  ++subjectIt;
		  }
		
		average += value * prob;
		weight += prob;
		}
	      }
	    }
	  
	  delete avgHistogram;
	  }
	    
	if ( weight > 0.0 ) 
	  {
	  filtered->Set( average / weight, offset );
	  } 
	else 
	  {
	  filtered->SetPaddingAt( offset );
	  }
	}
    }
  
  Progress::Done();
  
  return filtered;
}

} // namespace cmtk
