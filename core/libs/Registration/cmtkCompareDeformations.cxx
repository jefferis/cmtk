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

#include <cmtkCompareDeformations.h>

#include <cmtkVector3D.h>
#include <cmtkMathUtil.h>
#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

TypedArray* 
CompareDeformations::GetJacobianTTest
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
  Array<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  Array<Types::DataItem> valuesY( dataYsize );

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
	prob = MathUtil::TTest<Types::DataItem>( actualSizeX, valuesX, actualSizeY, valuesY, t, avgX, avgY );
	
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
CompareDeformations::GetPairedTTest
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
  Array<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  Array<Types::DataItem> valuesY( dataYsize );

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
      
      if ( actualSizeX == actualSizeY )
	{
	prob = MathUtil::PairedTTest<Types::DataItem>( actualSizeX, valuesX, valuesY, t, avgX, avgY );
	
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
CompareDeformations::GetPairedCorrelation
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
  Array<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  Array<Types::DataItem> valuesY( dataYsize );

  for ( unsigned int idx = 0; idx < length; ++idx ) 
    {
    correlationData->SetPaddingAt( idx );
    if ( pData )
      (*pData)->SetPaddingAt( idx );

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      unsigned int actualSize = 0;
      for ( unsigned int i = 0; i < dataXsize; ++i )
	if ( dataX[i]->Get( valuesX[actualSize], idx ) && dataY[i]->Get( valuesY[actualSize], idx ) )
	  ++actualSize;
      
      if ( actualSize )
	{
	Types::DataItem corr = MathUtil::Correlation<Types::DataItem>( actualSize, valuesX, valuesY );
	correlationData->Set(  corr, idx );
	if ( pData ) (*pData)->Set( MathUtil::ProbabilityFromTStat( MathUtil::TStatFromCorrelation( corr, actualSize-2 ), actualSize-2 ), idx );
	}
      }
    }
  
  return correlationData;
}

TypedArray* 
CompareDeformations::GetOneSampleTTest
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
  Array<Types::DataItem> valuesX( dataXsize );

  Types::DataItem t = 0, prob = 0, avgX = 0;
  for ( unsigned int idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) {
      unsigned int actualSizeX = 0;
      for ( unsigned int i = 0; i < dataXsize; ++i )
	if ( dataX[i]->Get( valuesX[actualSizeX], idx ) ) ++actualSizeX;
      
      if ( actualSizeX )
	{
	prob = MathUtil::TTest<Types::DataItem>( actualSizeX, valuesX, t, avgX );

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
CompareDeformations::GetHeritability
( std::vector<TypedArray::SmartPtr>& dataX, 
  std::vector<TypedArray::SmartPtr>& dataY,
  const TypedArray* mask )
{
  const size_t length = dataX[0]->GetDataSize();

  TypedArray* outData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataXsize = dataX.size();
  Array<float> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  Array<float> valuesY( dataYsize );

  for ( size_t idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) {
    }
  }
  
  return outData;
}

TypedArray* 
CompareDeformations::GetZScores
( std::vector<TypedArray::SmartPtr>& dataX,
  std::vector<TypedArray::SmartPtr>& dataY,
  const TypedArray* mask )
{
  const size_t length = dataX[0]->GetDataSize();

  TypedArray* outData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataXsize = dataX.size();
  Array<Types::DataItem> valuesX( dataXsize );
  const unsigned int dataYsize = dataY.size();
  Array<Types::DataItem> valuesY( dataYsize );

  Types::DataItem avgX, avgY, varX;

  for ( size_t idx = 0; idx < length; ++idx ) {

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
	avgX = MathUtil::Mean<Types::DataItem>( actualSizeX, valuesX );
	avgY = MathUtil::Mean<Types::DataItem>( actualSizeY, valuesY );
	
        varX = MathUtil::Variance<Types::DataItem>( actualSizeX, valuesX, avgX );

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
CompareDeformations::GetGeneticCovariance
( std::vector<TypedArray::SmartPtr>& dataMZ, 
  std::vector<TypedArray::SmartPtr>& dataDZ,
  const TypedArray* mask )
{
  const size_t length = dataMZ[0]->GetDataSize();

  TypedArray* outData = TypedArray::Create( TYPE_FLOAT, length );
  
  const unsigned int dataMZsize = dataMZ.size();
  Array<Types::DataItem> valuesMZ( dataMZsize );
  const unsigned int dataDZsize = dataDZ.size();
  Array<Types::DataItem> valuesDZ( dataDZsize );

  Types::DataItem avgMZ, avgDZ, varMZ, varDZ;

  for ( size_t idx = 0; idx < length; ++idx ) {

    Types::DataItem maskValue;
    if ( !mask || (mask->Get( maskValue, idx ) && (maskValue != 0)) ) 
      {
      unsigned int actualSizeMZ = 0;
      for ( unsigned int i = 0; i < dataMZsize; ++i )
	if ( dataMZ[i]->Get( valuesMZ[actualSizeMZ], idx ) ) ++actualSizeMZ;
      
      unsigned int actualSizeDZ = 0;
      for ( unsigned int i = 0; i < dataDZsize; ++i )
	if ( dataDZ[i]->Get( valuesDZ[actualSizeDZ], idx ) ) ++actualSizeDZ;

      if ( actualSizeMZ && actualSizeDZ )
	{
	avgMZ = MathUtil::Mean<Types::DataItem>( actualSizeMZ, valuesMZ );
	varMZ = MathUtil::Variance<Types::DataItem>( actualSizeMZ, valuesMZ, avgMZ );

	avgDZ = MathUtil::Mean<Types::DataItem>( actualSizeDZ, valuesDZ );
	varDZ = MathUtil::Variance<Types::DataItem>( actualSizeDZ, valuesDZ, avgDZ );

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
