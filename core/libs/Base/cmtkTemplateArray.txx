/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include <algorithm>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T>
double
TemplateArray<T>
::GetEntropy( const bool fractional, const int numberOfBins ) const 
{
  double entropy = 0;
  if ( fractional ) 
    {
    Histogram<double> histogram( numberOfBins );
    histogram.SetRange( this->GetRange() );
    
    for ( size_t idx = 0; idx < DataSize; ++idx )
      if ( !PaddingFlag || (Data[idx] != Padding) )
	histogram.IncrementFractional( histogram.ValueToBinFractional( Data[idx] ) );
    entropy = histogram.GetEntropy();
    } 
  else 
    {
    Histogram<unsigned int> histogram( numberOfBins );
    histogram.SetRange( this->GetRange() );
    
    for ( size_t idx = 0; idx < DataSize; ++idx )
      if ( !PaddingFlag || (Data[idx] != Padding) )
	histogram.Increment( histogram.ValueToBin( Data[idx] ) );
    entropy = histogram.GetEntropy();
    }
  return entropy;
}

template<class T> 
double
TemplateArray<T>
::GetEntropy( Histogram<unsigned int>& histogram ) const 
{
  histogram.Reset();
  for ( size_t idx = 0; idx < DataSize; ++idx )
    if ( !PaddingFlag || (Data[idx] != Padding) )
      histogram.Increment( histogram.ValueToBin( Data[idx] ) );
  return histogram.GetEntropy();
}

template<class T> double
TemplateArray<T>
::GetEntropy( Histogram<double>& histogram, const bool fractional ) const 
{
  histogram.Reset();
  if ( fractional ) 
    {
    for ( size_t idx = 0; idx < DataSize; ++idx )
      if ( !PaddingFlag || (Data[idx] != Padding) )
	histogram.IncrementFractional( histogram.ValueToBinFractional( Data[idx] ) );
    } 
  else
    {
    for ( size_t idx = 0; idx < DataSize; ++idx )
      if ( !PaddingFlag || (Data[idx] != Padding) )
	histogram.Increment( histogram.ValueToBin( Data[idx] ) );
    }
  return histogram.GetEntropy();
}

template<class T> double
TemplateArray<T>
::GetEntropy( Histogram<double>& histogram, const double* kernel, const size_t kernelRadius  ) const 
{
  histogram.Reset();
  for ( size_t idx = 0; idx < DataSize; ++idx )
    if ( !PaddingFlag || (Data[idx] != Padding) )
      histogram.AddWeightedSymmetricKernelFractional( histogram.ValueToBinFractional( Data[idx] ), kernelRadius, kernel );
  return histogram.GetEntropy();
}

template<class T>
const Types::DataItemRange
TemplateArray<T>
::GetRange() const 
{
  return Types::DataItemRange( this->GetRangeTemplate() );
}

template<class T> 
const Types::Range<T>
TemplateArray<T>
::GetRangeTemplate() 
  const 
{
  Types::Range<T> range( 0, 0 );
  
  // find first finite and non-Padding element
  size_t idx = 0;
  if ( PaddingFlag )
    {
    while ( (idx < DataSize) && ((Data[idx] == Padding) || !finite(Data[idx])) ) 
      ++idx;
    }
  else
    {
    while ( (idx < DataSize) && !finite(Data[idx]) ) 
      ++idx;
    }
  
  // didn't find any? return with error flag.
  if ( idx < DataSize) 
    {
    // otherwise: search for min/max from here
    range.m_LowerBound = range.m_UpperBound = Data[idx];

    if ( PaddingFlag ) 
      {
      for ( ; idx < DataSize; ++idx ) 
	{
	if ( (Data[idx] != Padding) && finite(Data[idx]) ) 
	  {
	  if (Data[idx] > range.m_UpperBound) 
	    range.m_UpperBound = Data[idx];
	  if (Data[idx] < range.m_LowerBound) 
	    range.m_LowerBound = Data[idx];
	  }
	}
      } 
    else
      {
      for ( ; idx < DataSize; ++idx ) 
	{
	if ( finite(Data[idx]) )
	  {
	  if (Data[idx] > range.m_UpperBound) 
	    range.m_UpperBound = Data[idx];
	  if (Data[idx] < range.m_LowerBound) 
	    range.m_LowerBound = Data[idx];
	  }
	}
      }
    }

  return range;
}

/**\todo This should somehow preserve PaddingFlag and Padding of the original 
 * array.
 */
template<class T>
TypedArray::SmartPtr
TemplateArray<T>
::Convert( const ScalarDataType dtype ) const 
{
  void* data = this->ConvertArray( dtype );
  return TypedArray::Create( dtype, data, DataSize, true /* freeArray */ );
}

template<class T>
void*
TemplateArray<T>
::ConvertSubArray 
( const ScalarDataType dtype, const size_t fromIdx, const size_t len ) const 
{
  void* data = Memory::AllocateArray<char>( len * TypeItemSize( dtype ) );
  this->ConvertSubArray( data, dtype, fromIdx, len );
  return data;
}

template<class T>
void
TemplateArray<T>
::GammaCorrection( const Types::DataItem gamma )
{
  if ( gamma > 0 ) 
    {
    Types::Range<T> range = this->GetRangeTemplate();
    
    const T diff = range.m_UpperBound - range.m_LowerBound;
    const double scale = 1.0 / diff;
    
#pragma omp parallel for if (DataSize>1e5)
    for ( size_t i = 0; i < DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) ) 
	{
	if ( Data[i] > range.m_LowerBound ) 
	  {
	  Data[i] = range.m_LowerBound + TypeTraits::Convert( diff * exp( log( scale * (Data[i]-range.m_LowerBound) ) / gamma ) );
	  }
	}
    }
}

template<class T>
void
TemplateArray<T>
::ApplyFunctionFloat( typename Self::FunctionTypeFloat f )
{
#pragma omp parallel for if (DataSize>1e5)
  for ( size_t i = 0; i < DataSize; ++i )
    if ( ! PaddingFlag || (Data[i] != Padding ) ) 
      {
      Data[i] = TypeTraits::Convert( f( (float)Data[i] ) );
      }
}

template<class T>
void
TemplateArray<T>
::ApplyFunctionDouble( typename Self::FunctionTypeDouble f )
{
#pragma omp parallel for if (DataSize>1e5)
  for ( size_t i = 0; i < DataSize; ++i )
    if ( ! PaddingFlag || (Data[i] != Padding ) ) 
      {
      Data[i] = TypeTraits::Convert( f( (double)Data[i] ) );
      }
}

template<class T>
void TemplateArray<T>
::ConvertSubArray 
( void *const destination, const ScalarDataType dtype, const size_t fromIdx, const size_t len ) const 
{
  if ( dtype == this->GetType() )
    memcpy( destination, Data + fromIdx, len * this->GetItemSize() );
  else 
    {
    switch ( dtype ) 
      {
      case TYPE_BYTE:
#pragma omp parallel for if (len>1e5)
	for ( size_t idx = 0; idx < len; ++idx )
	  ((byte*)destination)[idx] = DataTypeTraits<byte>::Convert( Data[ idx + fromIdx ] );
	break;
      case TYPE_CHAR:
#pragma omp parallel for if (len>1e5)
	for ( size_t idx = 0; idx < len; ++idx )
	  ((char*)destination)[idx] = DataTypeTraits<char>::Convert( Data[ idx + fromIdx ] );
	break;
      case TYPE_USHORT:
#pragma omp parallel for if (len>1e5)
	for ( size_t idx = 0; idx < len; ++idx )
	  ((unsigned short*)destination)[idx] = DataTypeTraits<unsigned short>::Convert( Data[ idx + fromIdx ] );
	break;
      case TYPE_SHORT:
#pragma omp parallel for if (len>1e5)
	for ( size_t idx = 0; idx < len; ++idx )
	  ((short*)destination)[idx] = DataTypeTraits<short>::Convert( Data[ idx + fromIdx ] );
	break;
      case TYPE_INT:
#pragma omp parallel for if (len>1e5)
	for ( size_t idx = 0; idx < len; ++idx )
	  ((int*)destination)[idx] = DataTypeTraits<int>::Convert( Data[ idx + fromIdx ] );
	break;
      case TYPE_FLOAT:
#pragma omp parallel for if (len>1e5)
	for ( size_t idx = 0; idx < len; ++idx )
	  ((float*)destination)[idx] = DataTypeTraits<float>::Convert( Data[ idx + fromIdx ] );
	break;
      case TYPE_DOUBLE:
#pragma omp parallel for if (len>1e5)
	for ( size_t idx = 0; idx < len; ++idx )
	  ((double*)destination)[idx] = DataTypeTraits<double>::Convert( Data[ idx + fromIdx ] );
	break;
      default:
	// unsupported / unknown data type. do nothing.
	break;
      }
    }
}

template<class T>
void
TemplateArray<T>
::ChangeEndianness() 
{
  size_t itemSize = this->GetItemSize();
  if ( itemSize < 2 ) return;
  size_t dataBytes = DataSize * itemSize;

  // f is the index of the FIRST byte of the current data item, l is the
  // index of that items LAST byte.
  size_t f, l;
  for ( f=0, l=itemSize-1; f<dataBytes; f+=itemSize,l+=itemSize )
    for ( size_t j=0; j<itemSize/2; ++j ) 
      {
      char d = ((char*)Data)[l-j];
      ((char*)Data)[l-j] = ((char*)Data)[f+j];
      ((char*)Data)[f+j] = d;
      }
}

template<class T>
void
TemplateArray<T>
::ApplyFunctionObject( const TypedArrayFunction& f )
{
#pragma omp parallel for
  for ( size_t i = 0; i < this->DataSize; ++i )
    {
    if ( !this->PaddingFlag || (this->Data[i] != this->Padding) )
      this->Data[i] = TypeTraits::Convert( f( this->Data[i] ) );
    }
}

template<class T>
void 
TemplateArray<T>
::BlockSet
( const Types::DataItem value, const size_t fromOffset, const size_t toOffset )
{
  T valueT = TypeTraits::Convert( value );

#pragma omp parallel for    
  for ( size_t i = fromOffset; i < toOffset; ++i )
    Data[i] = valueT;
}

} // namespace cmtk
