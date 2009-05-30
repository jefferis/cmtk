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

#include <cmtkScalarImageSimilarity.h>

#include <cmtkHistogram.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetMutualInformation
( const ScalarImage* image0, const ScalarImage* image1,
  ScalarImageSimilarityMemory *const memory )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  const TypedArray *data0 = image0->GetPixelData();
  const TypedArray *data1 = image1->GetPixelData();

  return TypedArraySimilarity::GetMutualInformation( data0, data1, memory );
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetNormalizedMutualInformation
( const ScalarImage* image0, const ScalarImage* image1, ScalarImageSimilarityMemory *const )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  const TypedArray *data0 = image0->GetPixelData();
  const TypedArray *data1 = image1->GetPixelData();

  return TypedArraySimilarity::GetNormalizedMutualInformation( data0, data1 );
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetMeanSquaredDifference
( const ScalarImage* image0, const ScalarImage* image1 )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  const TypedArray *data0 = image0->GetPixelData();
  const TypedArray *data1 = image1->GetPixelData();

  return TypedArraySimilarity::GetMeanSquaredDifference( data0, data1 );
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetCrossCorrelation
( const ScalarImage* image0, const ScalarImage* image1 )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  const TypedArray *data0 = image0->GetPixelData();
  const TypedArray *data1 = image1->GetPixelData();

  return TypedArraySimilarity::GetCrossCorrelation( data0, data1 );
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetGradientCorrelation
( const ScalarImage* image0, const ScalarImage* image1 )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  TypedArray::SmartPtr gradientX0( image0->GetSobelFiltered( CMTK_SCALARIMAGE_HORIZONTAL ) );
  TypedArray::SmartPtr gradientX1( image1->GetSobelFiltered( CMTK_SCALARIMAGE_HORIZONTAL ) );

  TypedArray::SmartPtr gradientY0( image0->GetSobelFiltered( CMTK_SCALARIMAGE_VERTICAL ) );
  TypedArray::SmartPtr gradientY1( image1->GetSobelFiltered( CMTK_SCALARIMAGE_VERTICAL ) );

  return TypedArraySimilarity::GetCrossCorrelation( gradientX0, gradientX1 ) + TypedArraySimilarity::GetCrossCorrelation( gradientY0, gradientY1 );
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetGradientDifference
( const ScalarImage* image0, const ScalarImage* image1, const ScalarImageSimilarity::ReturnType Ax, const ScalarImageSimilarity::ReturnType Ay )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  TypedArray::SmartPtr gradientX0( image0->GetSobelFiltered( CMTK_SCALARIMAGE_HORIZONTAL ) );
  TypedArray::SmartPtr gradientX1( image1->GetSobelFiltered( CMTK_SCALARIMAGE_HORIZONTAL ) );

  TypedArray::SmartPtr gradientY0( image0->GetSobelFiltered( CMTK_SCALARIMAGE_VERTICAL ) );
  TypedArray::SmartPtr gradientY1( image1->GetSobelFiltered( CMTK_SCALARIMAGE_VERTICAL ) );

  // absAx and absAy are the relative coefficients Ax and Ay multiplied with
  // the first gradient image's variances as described by Penney et al.
  Types::DataItem mean, var;
  gradientX0->GetStatistics( mean, var );
  ScalarImageSimilarity::ReturnType absAx = Ax * var;
  gradientY0->GetStatistics( mean, var );
  ScalarImageSimilarity::ReturnType absAy = Ay * var;

  Types::DataItem scaleFactor = 0;
  TypedArray::SmartPtr diffX( TypedArraySimilarity::GetDifferenceArray( gradientX0, gradientX1, scaleFactor ) );
  TypedArray::SmartPtr diffY( TypedArraySimilarity::GetDifferenceArray( gradientY0, gradientY1, scaleFactor ) );
  ScalarImageSimilarity::ReturnType resultX = 0, resultY = 0;
  
  unsigned int numberOfPixels = image0->GetNumberOfPixels();
  for ( unsigned int offset = 0; offset < numberOfPixels; ++offset ) 
    {
    Types::DataItem value;
    diffX->Get( value, offset );
    resultX += 1.0 / ( absAx + MathUtil::Square( value ) );
    diffY->Get( value, offset );
    resultY += 1.0 / ( absAy + MathUtil::Square( value ) );
    }
  
  return absAx * resultX + absAy * resultY;
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetPatternIntensity
( const ScalarImage* image0, const ScalarImage* image1, const ScalarImageSimilarity::ReturnType sigma, const unsigned int radius )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  ScalarImageSimilarity::ReturnType PI = 0;

  const int r2 = MathUtil::Square( radius );
  const ScalarImageSimilarity::ReturnType sigma2 = MathUtil::Square( sigma );

  // precompute desired offsets
  static std::vector<int> rows; 
  static std::vector<int> cols;

  static unsigned int lastRadius = 0;
  
  if ( radius != lastRadius )
    {
    lastRadius = radius;

    rows.clear();
    cols.clear();

    for ( int jj = -radius; jj <= (int)radius; ++jj ) 
      {
      for ( int ii = -radius; ii <= (int)radius; ++ii ) 
	{
	// are we inside circle?
	if ( ( MathUtil::Square(ii) + MathUtil::Square(jj) ) <= r2 ) 
	  {
	  rows.push_back(jj);
	  cols.push_back(ii);
	  }
	}
      }
    }
  
  // compute difference image
  Types::DataItem scaleFactor = 0;
  TypedArray::SmartPtr diff( TypedArraySimilarity::GetDifferenceArray( image0->GetPixelData(), image1->GetPixelData(), scaleFactor ) );
  
  // first, loop over all pixels in the images.
  const int dimsX = image0->GetDims(AXIS_X);
  const int dimsY = image0->GetDims(AXIS_Y);

  size_t offset = 0;
  for ( int j = 0; j < dimsY; ++j ) 
    {
    for ( int i = 0; i < dimsX; ++i, ++offset ) 
      {
      // get difference pixel in center of VOI around current pixel
      Types::DataItem centerValue;
      if ( diff->Get( centerValue, offset ) )
	{
	
	// loop over all pixels in square that encloses circular VOI
	for ( size_t rc=0; rc<rows.size(); rc++ ) 
	  {
	  const int jj = rows[rc];
	  const int ii = cols[rc];
	  
	  Types::DataItem VOIvalue;
	  // this is a hack really, since we're reading 2-D indexed data from
	  // a linear typed array. Because we generated this array as the
	  // subtraction of two 2-D images, however, this should work fine.
	  // doing it this way saves us the implementation of another,
	  // image-based rather than array-based subtraction function.
	  if ( (i+ii)>=0 && (i+ii)<dimsX && (j+jj)>=0 && (j+jj)<dimsY )
	    {
	    diff->Get( VOIvalue, offset + ii + dimsX * jj );
	    
	    // then update Pattern Intensity sum
	    PI += 1.0 / ( sigma2 + MathUtil::Square( centerValue - VOIvalue ) );
	    }
	  }
	}
      }
    }
  
  // normalize with sigma^2 that is common for all addends.
  PI *= sigma2;
  
  return PI;
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetDifferenceImageEntropy
( const ScalarImage* image0, const ScalarImage* image1 )
{
  Types::DataItem dummy;
  return GetDifferenceImageEntropy( image0, image1, dummy );
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetDifferenceImageEntropy
( const ScalarImage* image0, const ScalarImage* image1, Types::DataItem &scaleFactor )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;

  const TypedArray *data0 = image0->GetPixelData();
  const TypedArray *data1 = image1->GetPixelData();

  return TypedArraySimilarity::GetDifferenceArrayEntropy( data0, data1, scaleFactor );
}

ScalarImageSimilarity::ReturnType 
ScalarImageSimilarity::GetCorrelationRatio
( const ScalarImage* image0, const ScalarImage* image1 )
{
  if ( ! CheckImageDimensions( image0, image1 ) ) return CMTK_FLOAT_NAN;
  
  const TypedArray *data0 = image0->GetPixelData();
  const TypedArray *data1 = image1->GetPixelData();

  return TypedArraySimilarity::GetCorrelationRatio( data0, data1 );
}

bool 
ScalarImageSimilarity::CheckImageDimensions
( const ScalarImage* image0, const ScalarImage* image1 )
{
  if ( !image0 || !image1 ) return false;

  if ( !image0->GetPixelData() || !image1->GetPixelData() ) return false;

  return ( ( image0->GetDims()[0] == image1->GetDims()[0] ) && ( image0->GetDims()[1] == image1->GetDims()[1] ) );
}

} // namespace cmtk
