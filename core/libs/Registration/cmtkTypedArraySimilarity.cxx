/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkTypedArraySimilarity.h"

#include <Base/cmtkHistogram.h>
#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

TypedArraySimilarity::ReturnType 
TypedArraySimilarity::GetMutualInformation
( const TypedArray* array0, const TypedArray* array1,
  TypedArraySimilarityMemory *const memory )
{
  if ( ! CheckArrayDimensions( array0, array1 ) ) 
    return MathUtil::GetFloatNaN();

  size_t dataSize = array0->GetDataSize();

  JointHistogram<unsigned int>::SmartPtr histogram;
  if ( memory ) 
    {
    histogram = JointHistogram<unsigned int>::SmartPtr( memory->CreateHistogram( array0, array1 ) );
    }
  else
    {
    size_t numBins = std::max<unsigned>( std::min<unsigned>( static_cast<unsigned>( sqrt( (float)dataSize ) ), 128 ), 8 );
    
    histogram = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>( numBins, numBins ) );
    
    histogram->SetRangeX( array0->GetRange() );
    histogram->SetRangeY( array1->GetRange() );
    }

  Types::DataItem value0, value1;
  for ( unsigned int idx = 0; idx < dataSize; ++idx ) 
    {
    if ( array0->Get( value0, idx ) && array1->Get( value1, idx ) ) 
      {
      histogram->Increment( histogram->ValueToBinX( value0 ), histogram->ValueToBinY( value1 ) );
      }
    }
  
  return static_cast<TypedArraySimilarity::ReturnType>( histogram->GetMutualInformation( false ) );
}

TypedArraySimilarity::ReturnType 
TypedArraySimilarity::GetCorrelationRatio
 ( const TypedArray* array0, const TypedArray* array1 )
{
  // check if both images have same number of pixels.
  if ( ! CheckArrayDimensions( array0, array1 ) ) return MathUtil::GetFloatNaN();

  // determine reference image value range.
  const Types::DataItemRange range = array0->GetRange();

  // get pixel count and determine histogram size.
  const unsigned int dataSize = array0->GetDataSize();
  unsigned int numBins = std::max<unsigned>( std::min<unsigned>( static_cast<unsigned>( sqrt( (float)dataSize ) ), 128 ), 8 );

  // make sure number of classes doesn't exceed number of distinct values for
  // discrete data types.
  if ( (array0->GetType() != TYPE_FLOAT) && (array0->GetType() != TYPE_DOUBLE) ) 
    {
    numBins = std::min( numBins, static_cast<unsigned int>(range.Width()+1) );
    }
  
  // create histogram to count floating pixels in each reference class
  Histogram<unsigned int> histogram( numBins );

  // set value range for histogram to range of reference image.
  histogram.SetRange( range );

  // initialize arrays that hold the sums of all floating values and their
  // squares, separated by histogram classes of the reference image.
  double* sumJ = Memory::AllocateArray<double>( numBins );
  memset( sumJ, 0, numBins * sizeof( sumJ[0] ) );
  double* sumSquareJ = Memory::AllocateArray<double>( numBins );
  memset( sumSquareJ, 0, numBins * sizeof( sumSquareJ[0] ) );

  // sort all image intensities into data structures.
  Types::DataItem value0, value1;
  for ( unsigned int idx = 0; idx < dataSize; ++idx ) 
    {
    // for all valid voxel pairs
    if ( array0->Get( value0, idx ) && array1->Get( value1, idx ) ) 
      {
      // what's the reference histogram bin?
      unsigned int bin = histogram.ValueToBin( value0 );
      // count this sample
      histogram.Increment( bin );
      // add floating value to sum of values for this class
      sumJ[bin] += value1;
      // add squared floating value to sum of squared values for this class
      sumSquareJ[bin] += MathUtil::Square( value1 );
      }
    }
  
  double invSampleCount = 1.0 / histogram.SampleCount();
  // initialize variable for the weighted sum of the sigma^2 values over all
  // reference intensity classes.
  double sumSigmaSquare = 0;
  // run over all bins, i.e., reference classes
  for ( unsigned int j = 0; j < numBins; ++j ) 
    {
    // are there any values in the current class?
    if ( histogram[j] ) 
      {
      // compute mean floating value for this reference class
      double mu = sumJ[j] / histogram[j];
      // compute variance of floating values for this reference class
      double sigmaSq = ( mu*mu*histogram[j] - 2.0*mu*sumJ[j] + sumSquareJ[j] ) / histogram[j]; 
      // update sum over all classes with weighted sigma^2 for this class.
      sumSigmaSquare += (invSampleCount * histogram[j]) * sigmaSq;
      }
    }
  
  // get variance of complete floating image for normalization
  Types::DataItem sigmaSqJ, muJ;
  array1->GetStatistics( muJ, sigmaSqJ );

  Memory::DeleteArray( sumJ );
  Memory::DeleteArray( sumSquareJ );

  // return (supposedly) correlation ratio
  return 1.0 - (1.0 /  sigmaSqJ ) * sumSigmaSquare;
}

TypedArraySimilarity::ReturnType 
TypedArraySimilarity::GetNormalizedMutualInformation
( const TypedArray* array0, const TypedArray* array1,
  TypedArraySimilarityMemory *const memory )
{
  if ( ! CheckArrayDimensions( array0, array1 ) ) return MathUtil::GetFloatNaN();

  size_t dataSize = array0->GetDataSize();

  JointHistogram<unsigned int>::SmartPtr histogram;
  if ( memory ) 
    histogram = JointHistogram<unsigned int>::SmartPtr( memory->CreateHistogram( array0, array1 ) );
  else 
    {
    size_t numBins = std::max<unsigned>( std::min<unsigned>( static_cast<unsigned>( sqrt( (float)dataSize ) ), 128 ), 8 );
    
    histogram = JointHistogram<unsigned int>::SmartPtr( new JointHistogram<unsigned int>( numBins, numBins ) );
    histogram->SetRangeX( array0->GetRange() );
    histogram->SetRangeY( array1->GetRange() );
    }
  
  Types::DataItem value0, value1;
  for ( unsigned int idx = 0; idx < dataSize; ++idx ) 
    {
    if ( array0->Get( value0, idx ) && array1->Get( value1, idx ) ) 
      {
      histogram->Increment( histogram->ValueToBinX( value0 ), histogram->ValueToBinY( value1 ) );
      }
    }
  
  return static_cast<TypedArraySimilarity::ReturnType>( histogram->GetMutualInformation( true ) );
}

TypedArraySimilarity::ReturnType 
TypedArraySimilarity::GetMinusMeanSquaredDifference
( const TypedArray* array0, const TypedArray* array1 )
{
  if ( ! CheckArrayDimensions( array0, array1 ) ) return MathUtil::GetFloatNaN();

  unsigned int countPixels = 0;
  Types::DataItem pixel0, pixel1;
  Types::DataItem sumOfSquares = 0;

  unsigned int numberOfPixels = array0->GetDataSize();
  for ( unsigned int idx = 0; idx < numberOfPixels; ++idx ) 
    {
    if ( array0->Get( pixel0, idx ) && array1->Get( pixel1, idx ) ) 
      {
      sumOfSquares += MathUtil::Square( pixel0 - pixel1 );
      ++countPixels;
      }
    }
  
  if ( !countPixels )
    return MathUtil::GetFloatNaN();
  else
    return static_cast<TypedArraySimilarity::ReturnType>( -(sumOfSquares / (float)countPixels) );
}

TypedArraySimilarity::ReturnType
TypedArraySimilarity::GetPeakSignalToNoiseRatio
( const TypedArray* data, const TypedArray* signal )
{
  return -10.0 * log( -GetMinusMeanSquaredDifference( data, signal ) / signal->GetRange().Width() ) / log( 10.0 );
}

TypedArraySimilarity::ReturnType 
TypedArraySimilarity::GetCrossCorrelation
( const TypedArray* array0, const TypedArray* array1 )
{
  if ( ! CheckArrayDimensions( array0, array1 ) ) 
    return MathUtil::GetFloatNaN();

  Types::DataItem pixel0, pixel1;
  Types::DataItem sumOfProducts = 0, sumOfSquares0 = 0, sumOfSquares1 = 0;

  Types::DataItem mean0, mean1, dummy;
  array0->GetStatistics( mean0, dummy );
  array1->GetStatistics( mean1, dummy );

  unsigned int numberOfPixels = array0->GetDataSize();
  for ( unsigned int idx = 0; idx < numberOfPixels; ++idx ) 
    {
    if ( array0->Get( pixel0, idx ) && array1->Get( pixel1, idx ) ) 
      {
      sumOfProducts += (pixel0 - mean0) * (pixel1 - mean1);
      sumOfSquares0 += MathUtil::Square( pixel0 - mean0 );
      sumOfSquares1 += MathUtil::Square( pixel1 - mean1 );
      }
    }
  
  return sumOfProducts / ( sqrt( sumOfSquares0 ) * sqrt( sumOfSquares1 ) );
}

TypedArray::SmartPtr
TypedArraySimilarity::GetDifferenceArray
( const TypedArray* array0, const TypedArray* array1, Types::DataItem &scaleFactor )
{
  const size_t numberOfPixels = array0->GetDataSize();
  
  TypedArray::SmartPtr differenceArray = TypedArray::Create( GetSignedDataType( array0->GetType() ), numberOfPixels );
  
  Types::DataItem value0, value1;
  Types::DataItem ATA = 0.0, ATB = 0.0;
  for ( size_t i=0; i<numberOfPixels; i++) 
    {
    array0->Get( value0, i );
    ATA += (value0 * value0);
    
    array1->Get( value1, i );
    ATB += (value0 * value1);
    }

  // invert to get scale convention correct ( array0 = s*array1 )
  scaleFactor = ATA/ATB;

  Types::DataItem pixel0, pixel1;
  for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
    {
    if ( array0->Get( pixel0, idx ) && array1->Get( pixel1, idx ) ) 
      {
      differenceArray->Set( pixel0 - scaleFactor * pixel1, idx );
      }    
    }
  
  return differenceArray;
}

TypedArraySimilarity::ReturnType 
TypedArraySimilarity::GetDifferenceArrayEntropy
( const TypedArray* array0, const TypedArray* array1,
  Types::DataItem &scaleFactor )
{
  TypedArray::SmartPtr differenceArray( GetDifferenceArray( array0, array1, scaleFactor ) );

  return differenceArray->GetEntropy();
}

bool 
TypedArraySimilarity::CheckArrayDimensions
( const TypedArray* array0, const TypedArray* array1 )
{
  if ( !array0 || !array1 ) return false;

  return ( array0->GetDataSize() == array1->GetDataSize() );
}


TypedArraySimilarity::ReturnType
TypedArraySimilarity::GetOptimalScale
( const TypedArray* array0, const TypedArray* array1 )
{
  unsigned int dataSize = array0->GetDataSize();
  Types::DataItem value0, value1;

  TypedArraySimilarity::ReturnType ATA = 0.0;
  TypedArraySimilarity::ReturnType ATb = 0.0;

  for (unsigned int i=0; i<dataSize; i++) 
    {
    array0->Get( value0, i );
    ATA += (TypedArraySimilarity::ReturnType) (value0 * value0);
    
    array1->Get( value1, i );
    ATb += (TypedArraySimilarity::ReturnType) (value0 * value1);
    }
  
  return ATb/ATA;
}

} // namespace cmtk
