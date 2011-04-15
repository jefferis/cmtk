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
    
    for ( size_t idx = 0; idx < this->DataSize; ++idx )
      if ( !this->PaddingFlag || (this->Data[idx] != this->Padding) )
	histogram.IncrementFractional( histogram.ValueToBinFractional( this->Data[idx] ) );
    entropy = histogram.GetEntropy();
    } 
  else 
    {
    Histogram<unsigned int> histogram( numberOfBins );
    histogram.SetRange( this->GetRange() );
    
    for ( size_t idx = 0; idx < this->DataSize; ++idx )
      if ( !this->PaddingFlag || (this->Data[idx] != this->Padding) )
	histogram.Increment( histogram.ValueToBin( this->Data[idx] ) );
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
  for ( size_t idx = 0; idx < this->DataSize; ++idx )
    if ( !this->PaddingFlag || (this->Data[idx] != this->Padding) )
      histogram.Increment( histogram.ValueToBin( this->Data[idx] ) );
  return histogram.GetEntropy();
}

template<class T> double
TemplateArray<T>
::GetEntropy( Histogram<double>& histogram, const bool fractional ) const 
{
  histogram.Reset();
  if ( fractional ) 
    {
    for ( size_t idx = 0; idx < this->DataSize; ++idx )
      if ( !this->PaddingFlag || (this->Data[idx] != this->Padding) )
	histogram.IncrementFractional( histogram.ValueToBinFractional( this->Data[idx] ) );
    } 
  else
    {
    for ( size_t idx = 0; idx < this->DataSize; ++idx )
      if ( !this->PaddingFlag || (this->Data[idx] != this->Padding) )
	histogram.Increment( histogram.ValueToBin( this->Data[idx] ) );
    }
  return histogram.GetEntropy();
}

template<class T> double
TemplateArray<T>
::GetEntropy( Histogram<double>& histogram, const double* kernel, const size_t kernelRadius  ) const 
{
  histogram.Reset();
  for ( size_t idx = 0; idx < this->DataSize; ++idx )
    if ( !this->PaddingFlag || (this->Data[idx] != this->Padding) )
      histogram.AddWeightedSymmetricKernelFractional( histogram.ValueToBinFractional( this->Data[idx] ), kernelRadius, kernel );
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
  if ( this->PaddingFlag )
    {
    while ( (idx < this->DataSize) && ((this->Data[idx] == this->Padding) || !finite(this->Data[idx])) ) 
      ++idx;
    }
  else
    {
    while ( (idx < this->DataSize) && !finite(this->Data[idx]) ) 
      ++idx;
    }
  
  // didn't find any? return with error flag.
  if ( idx < this->DataSize) 
    {
    // otherwise: search for min/max from here
    range.m_LowerBound = range.m_UpperBound = this->Data[idx];

    if ( this->PaddingFlag ) 
      {
      for ( ; idx < this->DataSize; ++idx ) 
	{
	if ( (this->Data[idx] != this->Padding) && finite(this->Data[idx]) ) 
	  {
	  if (this->Data[idx] > range.m_UpperBound) 
	    range.m_UpperBound = this->Data[idx];
	  if (this->Data[idx] < range.m_LowerBound) 
	    range.m_LowerBound = this->Data[idx];
	  }
	}
      } 
    else
      {
      for ( ; idx < this->DataSize; ++idx ) 
	{
	if ( finite(this->Data[idx]) )
	  {
	  if (this->Data[idx] > range.m_UpperBound) 
	    range.m_UpperBound = this->Data[idx];
	  if (this->Data[idx] < range.m_LowerBound) 
	    range.m_LowerBound = this->Data[idx];
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
  return TypedArray::Create( dtype, data, this->DataSize, true /* freeArray */ );
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
    
#pragma omp parallel for if (this->DataSize>1e5)
    for ( int i = 0; i < static_cast<int>( this->DataSize ); ++i )
      if ( ! this->PaddingFlag || (this->Data[i] != this->Padding ) ) 
	{
	if ( this->Data[i] > range.m_LowerBound ) 
	  {
	  this->Data[i] = range.m_LowerBound + TypeTraits::Convert( diff * exp( log( scale * (this->Data[i]-range.m_LowerBound) ) / gamma ) );
	  }
	}
    }
}

template<class T>
void
TemplateArray<T>
::ApplyFunctionFloat( typename Self::FunctionTypeFloat f )
{
#pragma omp parallel for if (this->DataSize>1e5)
  for ( int i = 0; i < static_cast<int>( this->DataSize ); ++i )
    if ( ! this->PaddingFlag || (this->Data[i] != this->Padding ) ) 
      {
      this->Data[i] = TypeTraits::Convert( f( (float)this->Data[i] ) );
      }
}

template<class T>
void
TemplateArray<T>
::ApplyFunctionDouble( typename Self::FunctionTypeDouble f )
{
#pragma omp parallel for if (this->DataSize>1e5)
  for ( int i = 0; i < static_cast<int>( this->DataSize ); ++i )
    if ( ! this->PaddingFlag || (this->Data[i] != this->Padding ) ) 
      {
      this->Data[i] = TypeTraits::Convert( f( (double)this->Data[i] ) );
      }
}

template<class T>
void TemplateArray<T>
::ConvertSubArray 
( void *const destination, const ScalarDataType dtype, const size_t fromIdx, const size_t len ) const 
{
  if ( dtype == this->GetType() )
    memcpy( destination, this->Data + fromIdx, len * this->GetItemSize() );
  else 
    {
    switch ( dtype ) 
      {
      case TYPE_BYTE:
#pragma omp parallel for if (len>1e5)
	for ( int idx = 0; idx < static_cast<int>( len ); ++idx )
	  ((byte*)destination)[idx] = DataTypeTraits<byte>::Convert( this->Data[ idx + fromIdx ] );
	break;
      case TYPE_CHAR:
#pragma omp parallel for if (len>1e5)
	for ( int idx = 0; idx < static_cast<int>( len ); ++idx )
	  ((char*)destination)[idx] = DataTypeTraits<char>::Convert( this->Data[ idx + fromIdx ] );
	break;
      case TYPE_USHORT:
#pragma omp parallel for if (len>1e5)
	for ( int idx = 0; idx < static_cast<int>( len ); ++idx )
	  ((unsigned short*)destination)[idx] = DataTypeTraits<unsigned short>::Convert( this->Data[ idx + fromIdx ] );
	break;
      case TYPE_SHORT:
#pragma omp parallel for if (len>1e5)
	for ( int idx = 0; idx < static_cast<int>( len ); ++idx )
	  ((short*)destination)[idx] = DataTypeTraits<short>::Convert( this->Data[ idx + fromIdx ] );
	break;
      case TYPE_INT:
#pragma omp parallel for if (len>1e5)
	for ( int idx = 0; idx < static_cast<int>( len ); ++idx )
	  ((int*)destination)[idx] = DataTypeTraits<int>::Convert( this->Data[ idx + fromIdx ] );
	break;
      case TYPE_FLOAT:
#pragma omp parallel for if (len>1e5)
	for ( int idx = 0; idx < static_cast<int>( len ); ++idx )
	  ((float*)destination)[idx] = DataTypeTraits<float>::Convert( this->Data[ idx + fromIdx ] );
	break;
      case TYPE_DOUBLE:
#pragma omp parallel for if (len>1e5)
	for ( int idx = 0; idx < static_cast<int>( len ); ++idx )
	  ((double*)destination)[idx] = DataTypeTraits<double>::Convert( this->Data[ idx + fromIdx ] );
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
  size_t dataBytes = this->DataSize * itemSize;

  // f is the index of the FIRST byte of the current data item, l is the
  // index of that items LAST byte.
  size_t f, l;
  for ( f=0, l=itemSize-1; f<dataBytes; f+=itemSize,l+=itemSize )
    for ( size_t j=0; j<itemSize/2; ++j ) 
      {
      char d = ((char*)this->Data)[l-j];
      ((char*)this->Data)[l-j] = ((char*)this->Data)[f+j];
      ((char*)this->Data)[f+j] = d;
      }
}

template<class T>
void
TemplateArray<T>
::ApplyFunctionObject( const TypedArrayFunction& f )
{
#pragma omp parallel for
  for ( int i = 0; i < static_cast<int>( this->DataSize ); ++i )
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
  for ( int i = fromOffset; i < static_cast<int>( toOffset ); ++i )
    this->Data[i] = valueT;
}

} // namespace cmtk
