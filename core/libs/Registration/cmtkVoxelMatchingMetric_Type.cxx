/*
//
//  Copyright 1997-2012 Torsten Rohlfing
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

#include "cmtkVoxelMatchingMetric_Type.h"

#include <Base/cmtkJointHistogram.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class T,ScalarDataType DT>
void
VoxelMatchingMetric_Type<T,DT>::ImageData::Init
( const UniformVolume* volume )
{
  const TypedArray *srcArray = volume->GetData();
  DataArray = TypedArray::SmartPtr( srcArray->Convert( DT ) );
  Data = static_cast<T*>( DataArray->GetDataPtr() );
  NumberOfSamples = this->DataArray->GetDataSize();

  this->m_ValueRange = DataArray->GetRange();
  BinOffset = this->m_ValueRange.m_LowerBound;
  BinWidth = 1;

  if ( srcArray->GetPaddingFlag() )
    {
    Padding = DataTypeTraits<T>::Convert( srcArray->GetPaddingValue() );
    }
  else
    {
    Padding = DataTypeTraits<T>::ChoosePaddingValue();
    }
}

template<class T,ScalarDataType DT>
size_t
VoxelMatchingMetric_Type<T,DT>::ImageData::Init
( const UniformVolume* volume, const size_t defNumBins, const Types::DataItemRange& bounds )
{
  const TypedArray* data = volume->GetData();
  this->AllocDataArray( data );
  
  Types::DataItem value = 0, minValue = FLT_MAX, maxValue = -FLT_MAX;

  const DataGrid::IndexType cropFrom = volume->CropRegion().From();
  const DataGrid::IndexType cropTo = volume->CropRegion().To();
  const DataGrid::IndexType increments = volume->GetCropRegionIncrements();

  int offset = increments[0];
  for ( int z = cropFrom[2]; z < cropTo[2]; ++z, offset += increments[2] ) 
    {
    for ( int y = cropFrom[1]; y < cropTo[1]; ++y, offset += increments[1] ) 
      {
      for ( int x = cropFrom[0]; x < cropTo[0]; ++x, ++offset ) 
	{
	if ( data->Get( value, offset ) ) 
	  {
	  if ( value > maxValue ) maxValue = value;
	  if ( value < minValue ) minValue = value;
	  }
	}
      }
    }
  
  minValue = std::max( minValue, bounds.m_LowerBound );
  maxValue = std::min( maxValue, bounds.m_UpperBound );
  
  unsigned int defaultValue = 0;
  unsigned int numBins = defNumBins;

  if ( numBins != CMTK_HISTOGRAM_AUTOBINS ) 
    {
    BinOffset = minValue;
    BinWidth = ( maxValue - minValue ) / (numBins-1);
    double factor = 1.0 / BinWidth;
    
    for ( size_t idx = 0; idx < NumberOfSamples; ++idx ) 
      {
      if ( data->Get( value, idx ) ) 
	{
	value = std::max( std::min( value, maxValue ), minValue );
	Data[idx] = static_cast<T>( floor(factor * (value-BinOffset)) );
	} 
      else 
	{
	// point to extra bins at the end of each row/column for NULL data.
	Data[idx] = defaultValue;
	}
      }
    } 
  else 
    {
    switch ( data->GetDataClass() ) 
      {
      case DATACLASS_LABEL: 
      {
      numBins = 1 + static_cast<unsigned int>(maxValue-minValue);
      if ( numBins > 254 ) 
	{
	fprintf( stderr, "Fatal error: Cannot handle more than 254 different labels.\n" );
	exit( 1 );
	}
      
      BinOffset = 0;
      BinWidth = 1;
      
      for ( size_t idx = 0; idx < NumberOfSamples; ++idx ) 
	{
	if ( data->Get( value, idx ) )
	  Data[idx] = static_cast<T>( value - minValue );
	else
	  // point to extra bins at the end of each row/column for NULL data.
	  Data[idx] = defaultValue;
	}
      }
      break;
      default: // Handle everything else as grey-level data.
      case DATACLASS_GREY: 
      {
      numBins = JointHistogramBase::CalcNumBins( volume );
      BinOffset = minValue;
      BinWidth = ( maxValue - minValue ) / (numBins-1);
      double factor = 1.0 / BinWidth;
      
      for ( size_t idx = 0; idx < NumberOfSamples; ++idx ) 
	{
	if ( data->Get( value, idx ) ) 
	  {
	  value = std::max( std::min( value, maxValue ), minValue );
	  Data[idx] = static_cast<T>( floor(factor*(value-BinOffset)) );
	  } 
	else 
	  {
	  // point to extra bins at the end of each row/column for NULL data.
	  Data[idx] = defaultValue;
	  }
	}
      }
      break;
      } 
    }
  
  this->m_ValueRange = Types::DataItemRange( 0, static_cast<Types::DataItem>( numBins - 1 ) );
  
  return (Padding = numBins);
}

template<class T,ScalarDataType DT>
void
VoxelMatchingMetric_Type<T,DT>::ImageData::AllocDataArray
( const TypedArray* templateArray )
{
  NumberOfSamples = templateArray->GetDataSize();
  DataArray = TypedArray::SmartPtr( TypedArray::Create( DT, NumberOfSamples ) );
  Data = static_cast<T*>( DataArray->GetDataPtr() );
}

template<class T,ScalarDataType DT>
void
VoxelMatchingMetric_Type<T,DT>::ImageData::PrecomputeIncrements
( const UniformVolume* volume )
{
  this->ImageDims = volume->GetDims();
  // pre-compute relative offsets for tri-linear model interpolation
  nextJ = volume->GetDims()[0];
  nextK = nextJ * volume->GetDims()[1];
  nextIJ = nextJ + 1;
  nextIK = nextK + 1;
  nextJK = nextK + nextJ;
  nextIJK = nextJK + 1;
}

// instantiate necessary templates.
template class VoxelMatchingMetric_Type<byte,TYPE_BYTE>;
template class VoxelMatchingMetric_Type<short,TYPE_SHORT>;

} // namespace cmtk
