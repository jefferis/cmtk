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

#include <cmtkHypothesisTests.h>

#include <cmtkVector3D.h>
#include <cmtkMathUtil.h>
#include <cmtkConsole.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

TypedArray* 
HypothesisTests::GetUnpairedTwoTailedTTest
( std::vector<TypedArray::SmartPtr>& dataX, 
  std::vector<TypedArray::SmartPtr>& dataY,
  TypedArray** tstatData, TypedArray** avgXData, 
  TypedArray** avgYData, const TypedArray* mask )
{
  const unsigned int length = dataX[0]->GetDataSize();

  TypedArray* probData = TypedArray::Create( TYPE_FLOAT, length );

  if ( tstatData )
    *tstatData = TypedArray::Create( TYPE_FLOAT, length );
  
  if ( avgXData )
    *avgXData = TypedArray::Create( TYPE_FLOAT, length );
  
  if ( avgYData )
    *avgYData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataXsize = dataX.size();
  std::vector<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  std::vector<Types::DataItem> valuesY( dataYsize );

  Types::DataItem t = 0, prob = 0, avgX = 0, avgY = 0;
  for ( unsigned int idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      unsigned int actualSizeX = 0;
      for ( unsigned int i = 0; i < dataXsize; ++i )
	if ( dataX[i]->Get( valuesX[actualSizeX], idx ) ) ++actualSizeX;
     
      unsigned int actualSizeY = 0;
      for ( unsigned int i = 0; i < dataYsize; ++i )
	if ( dataY[i]->Get( valuesY[actualSizeY], idx ) ) ++actualSizeY;
      
      if ( actualSizeX && actualSizeY )
	{
	prob = MathUtil::TTest<Types::DataItem>( valuesX, valuesY, t, avgX, avgY );
	
	if ( (prob < 0) || (prob>1) )
	  {
	  fprintf( stderr, "t = %f\tp = %f\n", t, prob );
	  }
	prob = 1.0 - prob; // convert probability to significance
	}
      else
	{
	t = prob = 0;
	}
     
      if ( tstatData ) (*tstatData)->Set( t, idx );
      if ( avgXData ) (*avgXData)->Set( avgX, idx );
      if ( avgYData ) (*avgYData)->Set( avgY, idx );
      
      if ( avgX > avgY )
	probData->Set(  prob, idx );
      else
	probData->Set( -prob, idx );
      } 
    else 
      {
      probData->SetPaddingAt( idx );
      if ( tstatData ) (*tstatData)->SetPaddingAt( idx );
      if ( avgXData ) (*avgXData)->SetPaddingAt( idx );
      if ( avgYData ) (*avgYData)->SetPaddingAt( idx );
      }
  }
  
  return probData;
}

TypedArray* 
HypothesisTests::GetPairedTwoTailedTTest
( std::vector<TypedArray::SmartPtr>& dataX, 
  std::vector<TypedArray::SmartPtr>& dataY,
  TypedArray** tstatData, TypedArray** avgXData, 
  TypedArray** avgYData, const TypedArray* mask )
{
  if ( dataX.size() != dataY.size() )
    {
    StdErr << "Cannot perform paired t-test if numbers of X and Y samples isn't equal.\n";
    return NULL;
    }
      
  const unsigned int length = dataX[0]->GetDataSize();

  TypedArray* probData = TypedArray::Create( TYPE_FLOAT, length );

  if ( tstatData )
    *tstatData = TypedArray::Create( TYPE_FLOAT, length );
  
  if ( avgXData )
    *avgXData = TypedArray::Create( TYPE_FLOAT, length );
  
  if ( avgYData )
    *avgYData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataXsize = dataX.size();
  std::vector<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  std::vector<Types::DataItem> valuesY( dataYsize );

  Types::DataItem t = 0, prob = 0, avgX = 0, avgY = 0;
  for ( unsigned int idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      valuesX.resize( dataXsize );
      unsigned int actualSizeX = 0;
      for ( unsigned int i = 0; i < dataXsize; ++i )
	if ( dataX[i]->Get( valuesX[actualSizeX], idx ) ) ++actualSizeX;
      
      valuesY.resize( dataYsize );
      unsigned int actualSizeY = 0;
      for ( unsigned int i = 0; i < dataYsize; ++i )
	if ( dataY[i]->Get( valuesY[actualSizeY], idx ) ) ++actualSizeY;
      
      if ( actualSizeX == actualSizeY )
	{
	valuesX.resize( actualSizeX );
	valuesY.resize( actualSizeY );

	prob = MathUtil::PairedTTest<Types::DataItem>( valuesX, valuesY, t, avgX, avgY );
	
	if ( (prob < 0) || (prob>1) )
	  {
	  fprintf( stderr, "t = %f\tp = %f\n", t, prob );
	  }
	prob = 1.0 - prob; // convert probability to significance
	}
      else
	{
	t = prob = 0;
	}
      
      if ( tstatData ) (*tstatData)->Set( t, idx );
      if ( avgXData ) (*avgXData)->Set( avgX, idx );
      if ( avgYData ) (*avgYData)->Set( avgY, idx );
      
      if ( avgX > avgY )
	probData->Set(  prob, idx );
      else
	probData->Set( -prob, idx );
      } 
    else 
      {
      probData->SetPaddingAt( idx );
      if ( tstatData ) (*tstatData)->SetPaddingAt( idx );
      if ( avgXData ) (*avgXData)->SetPaddingAt( idx );
      if ( avgYData ) (*avgYData)->SetPaddingAt( idx );
      }
  }
  
  return probData;
}

TypedArray* 
HypothesisTests::GetPairedCorrelation
( std::vector<TypedArray::SmartPtr>& dataX, std::vector<TypedArray::SmartPtr>& dataY, TypedArray** pData, const TypedArray* mask )
{
  if ( dataX.size() != dataY.size() )
    {
    StdErr << "Cannot perform paired correlation if numbers of X and Y samples isn't equal.\n";
    return NULL;
    }
  
  const unsigned int length = dataX[0]->GetDataSize();

  TypedArray* correlationData = TypedArray::Create( TYPE_FLOAT, length );
  if ( pData )
    *pData = TypedArray::Create( TYPE_FLOAT, length );

  const unsigned int dataXsize = dataX.size();
  std::vector<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  std::vector<Types::DataItem> valuesY( dataYsize );

  for ( unsigned int idx = 0; idx < length; ++idx ) 
    {
    correlationData->SetPaddingAt( idx );
    if ( pData )
      (*pData)->SetPaddingAt( idx );

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      valuesX.resize( dataXsize );
      valuesY.resize( dataXsize );

      unsigned int actualSize = 0;
      for ( unsigned int i = 0; i < dataXsize; ++i )
	if ( dataX[i]->Get( valuesX[actualSize], idx ) && dataY[i]->Get( valuesY[actualSize], idx ) )
	  ++actualSize;
      
      if ( actualSize )
	{
	valuesX.resize( actualSize );
	valuesY.resize( actualSize );

	Types::DataItem corr = MathUtil::Correlation<Types::DataItem>( valuesX, valuesY );
	correlationData->Set(  corr, idx );
	if ( pData ) 
	  (*pData)->Set( MathUtil::ProbabilityFromTStat( MathUtil::TStatFromCorrelation( corr, actualSize-2 ), actualSize-2 ), idx );
	}
      }
    }
  
  return correlationData;
}

TypedArray* 
HypothesisTests::GetOneSampleTTest
( std::vector<TypedArray::SmartPtr>& dataX, 
  TypedArray** tstatData, TypedArray** avgXData, 
  const TypedArray* mask )
{
  const unsigned int length = dataX[0]->GetDataSize();

  TypedArray* probData = TypedArray::Create( TYPE_FLOAT, length );

  if ( tstatData )
    *tstatData = TypedArray::Create( TYPE_FLOAT, length );
  
  if ( avgXData )
    *avgXData = TypedArray::Create( TYPE_FLOAT, length );

  const unsigned int dataXsize = dataX.size();
  std::vector<Types::DataItem> valuesX( dataXsize );

  Types::DataItem t = 0, prob = 0, avgX = 0;
  for ( unsigned int idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      valuesX.resize( dataXsize );
      unsigned int actualSizeX = 0;
      for ( unsigned int i = 0; i < dataXsize; ++i )
	if ( dataX[i]->Get( valuesX[actualSizeX], idx ) ) 
	  ++actualSizeX;
      
      if ( actualSizeX )
	{
	valuesX.resize( actualSizeX );
	prob = MathUtil::TTest<Types::DataItem>( valuesX, t, avgX );

	if ( (prob < 0) || (prob>1) )
	  {
	  fprintf( stderr, "t = %f\tp = %f\n", t, prob );
	  }
	prob = 1.0 - prob; // convert probability to significance
	}
      else
	{
	t = prob = 0;
	}
      
      if ( tstatData ) (*tstatData)->Set( t, idx );
      if ( avgXData ) (*avgXData)->Set( avgX, idx );
      
      if ( avgX > 0 )
	probData->Set(  -prob, idx );
      else
	probData->Set( +prob, idx );
    } else {
      probData->SetPaddingAt( idx );
      if ( tstatData ) (*tstatData)->SetPaddingAt( idx );
      if ( avgXData ) (*avgXData)->SetPaddingAt( idx );
    }
  }
  
  return probData;
}

TypedArray* 
HypothesisTests::GetHeritability
( std::vector<TypedArray::SmartPtr>& dataX, 
  std::vector<TypedArray::SmartPtr>& dataY,
  const TypedArray* mask )
{
  const size_t length = dataX[0]->GetDataSize();

  TypedArray* outData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataXsize = dataX.size();
  std::vector<float> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  std::vector<float> valuesY( dataYsize );

  for ( size_t idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) {
    }
  }
  
  return outData;
}

TypedArray* 
HypothesisTests::GetZScores
( std::vector<TypedArray::SmartPtr>& dataX,
  std::vector<TypedArray::SmartPtr>& dataY,
  const TypedArray* mask )
{
  const size_t length = dataX[0]->GetDataSize();

  TypedArray* outData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataXsize = dataX.size();
  std::vector<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  std::vector<Types::DataItem> valuesY( dataYsize );

  Types::DataItem avgX, avgY, varX;

  for ( size_t idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      valuesX.resize( dataXsize );
      unsigned int actualSizeX = 0;
      for ( unsigned int i = 0; i < dataXsize; ++i )
	if ( dataX[i]->Get( valuesX[actualSizeX], idx ) ) ++actualSizeX;
      
      valuesY.resize( dataYsize );
      unsigned int actualSizeY = 0;
      for ( unsigned int i = 0; i < dataYsize; ++i )
	if ( dataY[i]->Get( valuesY[actualSizeY], idx ) ) ++actualSizeY;

      if ( actualSizeX && actualSizeY )
	{
	valuesX.resize( actualSizeX );
	avgX = MathUtil::Mean<Types::DataItem>( valuesX );
	valuesY.resize( actualSizeY );
	avgY = MathUtil::Mean<Types::DataItem>( valuesY );
	
        varX = MathUtil::Variance<Types::DataItem>( valuesX, avgX );

	outData->Set( (avgY - avgX) / sqrt( varX ), idx );
	}
      else
	{
	outData->SetPaddingAt( idx );
	}
      }
    else
      {
      outData->SetPaddingAt( idx );
      }
  }
  
  return outData;
}

TypedArray* 
HypothesisTests::GetGeneticCovariance
( std::vector<TypedArray::SmartPtr>& dataMZ, 
  std::vector<TypedArray::SmartPtr>& dataDZ,
  const TypedArray* mask )
{
  const size_t length = dataMZ[0]->GetDataSize();

  TypedArray* outData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataMZsize = dataMZ.size();
  std::vector<Types::DataItem> valuesMZ( dataMZsize );
  const unsigned int dataDZsize = dataDZ.size();
  std::vector<Types::DataItem> valuesDZ( dataDZsize );

  Types::DataItem avgMZ, avgDZ, varMZ, varDZ;

  for ( size_t idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      valuesMZ.resize( dataMZsize );
      unsigned int actualSizeMZ = 0;
      for ( unsigned int i = 0; i < dataMZsize; ++i )
	if ( dataMZ[i]->Get( valuesMZ[actualSizeMZ], idx ) ) ++actualSizeMZ;
      
      valuesDZ.resize( dataDZsize );
      unsigned int actualSizeDZ = 0;
      for ( unsigned int i = 0; i < dataDZsize; ++i )
	if ( dataDZ[i]->Get( valuesDZ[actualSizeDZ], idx ) ) ++actualSizeDZ;

      if ( actualSizeMZ && actualSizeDZ )
	{
	valuesMZ.resize( actualSizeMZ );
	avgMZ = MathUtil::Mean<Types::DataItem>( valuesMZ );
	varMZ = MathUtil::Variance<Types::DataItem>( valuesMZ, avgMZ );

	valuesDZ.resize( actualSizeDZ );
	avgDZ = MathUtil::Mean<Types::DataItem>( valuesDZ );
	varDZ = MathUtil::Variance<Types::DataItem>( valuesDZ, avgDZ );

	outData->Set( varMZ / avgMZ - varDZ / avgDZ, idx );
	}
      else
	{
	outData->SetPaddingAt( idx );
	}
      }
    else
      {
      outData->SetPaddingAt( idx );
      }
  }
  
  return outData;
}

} // namespace cmtk
