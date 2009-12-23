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

#include <algorithm>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

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
    
    Types::DataItem min, max;
    this->GetRange( min, max );
    histogram.SetRange( min, max );
    
    for ( size_t idx = 0; idx < DataSize; ++idx )
      if ( !PaddingFlag || (Data[idx] != Padding) )
	histogram.IncrementFractional( histogram.ValueToBinFractional( Data[idx] ) );
    entropy = histogram.GetEntropy();
    } 
  else 
    {
    Histogram<unsigned int> histogram( numberOfBins );
    
    Types::DataItem min, max;
    this->GetRange( min, max );
    histogram.SetRange( min, max );
    
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

template<class T> bool
TemplateArray<T>
::GetRange ( Types::DataItem& min, Types::DataItem& max ) const 
{
  T minT = 0, maxT = 0;
  if ( this->GetRangeTemplate( minT, maxT ) ) 
    {
    min = static_cast<Types::DataItem>( minT );
    max = static_cast<Types::DataItem>( maxT );
    return true;
    }
  return false;
}

template<class T> bool
TemplateArray<T>
::GetRangeTemplate( T& min, T& max ) 
  const 
{
  // find first finite and non-Padding element
  size_t idx = 0;
  if ( PaddingFlag )
    {
    while ( (idx < DataSize) && ((Data[idx] == Padding) || !finite(Data[idx])) ) 
      ++idx;
    }
  else
    {
    while ( !finite(Data[idx]) ) 
      ++idx;
    }

  // didn't find any? return with error flag.
  if ( !(idx < DataSize) ) 
    {
    return false;
    }

  // otherwise: search for min/max from here
  min = Data[idx];
  max = min;

  if ( PaddingFlag ) 
    {
    for ( ; idx < DataSize; ++idx ) 
      {
      if ( (Data[idx] != Padding) && finite(Data[idx]) ) 
	{
	if (Data[idx] > max) max = Data[idx];
	if (Data[idx] < min) min = Data[idx];
	}
      }
    } 
  else
    {
    for ( ; idx < DataSize; ++idx ) 
      {
      if ( finite(Data[idx]) )
	{
	if (Data[idx] > max) max = Data[idx];
	if (Data[idx] < min) min = Data[idx];
	}
      }
    }
  return true;
}

/**\todo This should somehow preserve PaddingFlag and Padding of the original 
 * array.
 */
template<class T>
TypedArray* 
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
  void* data = malloc( len * TypeItemSize( dtype ) );
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
    T min, max;
    this->GetRangeTemplate( min, max );
    
    T range = max - min;
    const double scale = 1.0 / range;

#pragma omp parallel for if (DataSize>1e5)
    for ( size_t i = 0; i < DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) ) 
	{
	if ( Data[i] > min ) 
	  {
	  Data[i] = min + TypeTraits::Convert( range * exp( log( scale * (Data[i]-min) ) / gamma ) );
	  }
	}
    }
}

template<class T>
void
TemplateArray<T>
::ApplyFunctionFloat( Self::FunctionTypeFloat f )
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
::ApplyFunctionDouble( Self::FunctionTypeDouble f )
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
::HistogramEqualization( const int numberOfLevels ) 
{
  std::vector<unsigned int> histogram( numberOfLevels );
  std::fill( histogram.begin(), histogram.end(), 0 );

  // find original value range
  T min = Data[0], max = Data[0];
  for ( size_t i = 1; i < DataSize; ++i )
    if ( !PaddingFlag || (Data[i] != Padding) ) 
      {
      if ( Data[i] > max ) max = Data[i];
      if ( Data[i] < min ) min = Data[i];
      }

  // generate histogram
  const double scale = (1.0 * (numberOfLevels-1)) / (max-min);
  for ( size_t i = 0; i < DataSize; ++i ) 
    {
    if ( !PaddingFlag || (Data[i] != Padding) ) 
      ++histogram[static_cast<unsigned int>(scale*(Data[i]-min))];
    }
  
  // transform into cumulative histogram
  histogram[0] = 0; // this effectively stretches the distribution
  for ( int l = 1; l < numberOfLevels; ++l ) 
    {
    histogram[l] += histogram[l-1];
    }
  
  // apply equalization tranformation
  const double stretchScale = 1.0 * (max-min) / histogram[numberOfLevels-1];
#pragma omp parallel for
  for ( size_t i = 0; i < DataSize; ++i ) 
    {
    if ( !PaddingFlag || (Data[i] != Padding) ) 
      {
      const double equalized = stretchScale * histogram[static_cast<unsigned int>( scale * (Data[i]-min) )];
      Data[i] = min + std::max<T>( 0, static_cast<T>( this->ConvertItem( equalized ) ) );
      }
    }
}

template<class T>
void
TemplateArray<T>
::MatchHistogramToReference( const TypedArray* referenceArray, const unsigned int numberOfBins )
{
  Histogram<unsigned int>::SmartPtr referenceHistogram( referenceArray->GetHistogram( numberOfBins ) );
  referenceHistogram->ConvertToCumulative();
  
  std::vector<double> normalizedRefHistogram( numberOfBins );
  for ( size_t l = 0; l < numberOfBins; ++l )
    {
    normalizedRefHistogram[l] = 1.0 * referenceHistogram->GetBin(l) / referenceHistogram->GetBin(numberOfBins-1);
    }
  
  Histogram<unsigned int>::SmartPtr movingHistogram( this->GetHistogram( numberOfBins ) );
  movingHistogram->ConvertToCumulative();

  std::vector<double> normalizedMovHistogram( numberOfBins );
  for ( size_t l = 0; l < numberOfBins; ++l )
    {
    normalizedMovHistogram[l] =  1.0 * movingHistogram->GetBin(l) / movingHistogram->GetBin(numberOfBins-1);
    }
  
  std::vector<unsigned int> lookup( numberOfBins );
  size_t j = 0;
  for ( size_t i = 0; i < numberOfBins; ++i )
    {
    while ((j < numberOfBins) && (normalizedRefHistogram[j] < normalizedMovHistogram[i]))
      {
      ++j;
      }
    lookup[i] = j;
    }
  
  for ( size_t i = 0; i < this->DataSize; ++i )
    {
    if ( !this->PaddingFlag || (this->Data[i] != this->Padding) )
      this->Data[i] = static_cast<T>( referenceHistogram->BinToValue( lookup[ movingHistogram->ValueToBin( this->Data[i] ) ] ) );
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
