/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include <cmtkTypedArray.h>
#include <cmtkTemplateArray.h>

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
TypedArray
::RescaleToRange( const Types::DataItem min, const Types::DataItem max )
{
  Types::DataItem cMin, cMax;
  this->GetRange( cMin, cMax );

  const Types::DataItem scale = (max-min) / (cMax-cMin);
  const Types::DataItem offset = min - (cMin * scale);

  this->Rescale( scale, offset );
}

void
TypedArray
::BlockSwap
( const size_t fromOffset, const size_t toOffset, const size_t blockLength )
{
  char buffer[2048];

  size_t itemSize = this->GetItemSize();
  char *dataPtr = static_cast<char*>( this->GetDataPtr() );

  size_t bytesToGo = itemSize * blockLength;
  char *fromPtr = dataPtr + itemSize * fromOffset;
  char *toPtr = dataPtr + itemSize * toOffset;

  while ( bytesToGo > sizeof( buffer ) ) 
    {
    memcpy( buffer, toPtr, sizeof( buffer ) );
    memcpy( toPtr, fromPtr, sizeof( buffer ) );
    memcpy( fromPtr, buffer, sizeof( buffer ) );
    
    fromPtr += sizeof( buffer );
    toPtr += sizeof( buffer );
    bytesToGo -= sizeof( buffer );
    }

  if ( bytesToGo ) 
    {
    memcpy( buffer, toPtr, bytesToGo );
    memcpy( toPtr, fromPtr, bytesToGo );
    memcpy( fromPtr, buffer, bytesToGo );
    }
}

void
TypedArray
::BlockReverse
( const size_t fromOffset, const size_t blockLength )
{
  // we get into trouble here for item sizes > 128 bits.
  char buffer[16];

  size_t itemSize = this->GetItemSize();
  char *dataPtr = static_cast<char*>( this->GetDataPtr() );

  char *fromPtr = dataPtr + fromOffset * itemSize;
  char *toPtr = fromPtr + (blockLength-1) * itemSize;

  for ( size_t itemsToGo = blockLength / 2; itemsToGo; --itemsToGo ) 
    {
    memcpy( buffer, toPtr, itemSize );
    memcpy( toPtr, fromPtr, itemSize );
    memcpy( fromPtr, buffer, itemSize );
    
    fromPtr += itemSize;
    toPtr -= itemSize;
    }
}

TypedArray*
TypedArray
::Create
( const ScalarDataType dtype, void *const data, const size_t size,
  const bool freeArray, const bool paddingFlag, const void* paddingData ) 
{
  TypedArray *instance = NULL;
  switch (dtype) 
    {
    case TYPE_BYTE: 
      instance = new ByteArray( data, size, freeArray, paddingFlag, paddingData );
      break;
    case TYPE_CHAR: 
      instance = new CharArray( data, size, freeArray, paddingFlag, paddingData );
      break;
    case TYPE_SHORT: 
      instance = new ShortArray( data, size, freeArray, paddingFlag, paddingData );
      break;
    case TYPE_USHORT: 
      instance = new UShortArray( data, size, freeArray, paddingFlag, paddingData );
      break;
    case TYPE_INT: 
      instance = new IntArray( data, size, freeArray, paddingFlag, paddingData );
      break;
    case TYPE_FLOAT: 
      instance = new FloatArray( data, size, freeArray, paddingFlag, paddingData );
      break;
    case TYPE_DOUBLE: 
      instance = new DoubleArray( data, size, freeArray, paddingFlag, paddingData );
      break;
    default:
      fprintf(stderr,"TypedArray::Create - Data type %d unknown.",dtype);
      abort();
    }
  
  return instance;
}

TypedArray*
TypedArray
::Create( const ScalarDataType dtype, const size_t size )
{
  TypedArray *instance = NULL;
  switch ( dtype ) 
    {
    case TYPE_BYTE:   
      instance = new ByteArray( size );
      break;
    case TYPE_CHAR:   
      instance = new CharArray( size );
      break;
    case TYPE_SHORT:  
      instance = new ShortArray( size );
      break;
    case TYPE_USHORT:  
      instance = new UShortArray( size );
      break;
    case TYPE_INT:    
      instance = new IntArray( size );
      break;
    case TYPE_FLOAT:  
      instance = new FloatArray( size );
      break;
    case TYPE_DOUBLE:
      instance = new DoubleArray( size );
      break;
    default:
      fprintf( stderr, "TypedArray::Create - Data type %d unknown.", dtype );
      abort();
    }
  
  return instance;
}

void 
TypedArray
::PruneHistogram
( const bool pruneHi, const bool pruneLo, const size_t numberOfBinsTarget, const size_t numberOfBinsInternal )
{
  Histogram<unsigned int>::SmartPtr originalHistogram( this->GetHistogram( numberOfBinsInternal ) );
  
  const size_t oneBinFraction = this->GetDataSize() / numberOfBinsTarget;
  size_t accumulatedNumberOfSamples = 0;

  Types::DataItem min, max;
  this->GetRange( min, max );
  const Types::DataItem originalMax = max;

  if ( pruneHi )
    {
    for ( size_t binIdx = numberOfBinsInternal-1; binIdx > 0; --binIdx )
      {
      accumulatedNumberOfSamples += (*originalHistogram)[binIdx];
      if ( accumulatedNumberOfSamples > oneBinFraction )
	{
	max = min + (max-min)/1024*binIdx;
	break;
	}
      }
    }

  if ( pruneLo )
    {
    for ( size_t binIdx = 0; binIdx < numberOfBinsInternal; ++binIdx )
      {
      accumulatedNumberOfSamples += (*originalHistogram)[binIdx];
      if ( accumulatedNumberOfSamples > oneBinFraction )
	{
	min = originalMax - (originalMax-min)/numberOfBinsInternal*binIdx;
	break;
	}
      }
    }
  
  this->Threshold( min, max );
}

} // namespace cmtk
