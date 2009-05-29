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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkTypedArrayRLE.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class T>
TemplateArrayRLE<T>::TemplateArrayRLE
( const TypedArray* dataArray, const unsigned cacheBlockLength ) 
{
  if ( cacheBlockLength )
    CacheBlockLength = cacheBlockLength;
  else
    CacheBlockLength = Self::DefaultCacheBlockLength;
  
  if ( dataArray->GetType() != this->GetType() )
    this->Data = static_cast<T*>( dataArray->ConvertArray( this->GetType() ) );
  else
    this->Data = // this is okay since we really only use Data for reading.
      const_cast<T*>( reinterpret_cast<const T*>( dataArray->GetDataPtr() ) );
  
  // get length of RLE data and create data array
  this->DataSize = dataArray->GetDataSize();
  NumberEntriesRLE = this->RequiredLengthRLE();
  DataRLE = Memory::AllocateArray<RLEentry>( NumberEntriesRLE+1 ); // allocate +1 for guard element
  
  // fill RLE data from original data
  typename Self::RLEentry *entry = DataRLE;
  entry->Value = this->Data[0];
  entry->Multiplicity = 1;
  
  // also create access cache along the way
  AccessCache = // alloc. +1 for incomplete block and +1 for guard element
    Memory::AllocateArray<typename Self::AccessCacheEntry>(  this->DataSize / CacheBlockLength + 2 );
  typename Self::AccessCacheEntry *cache = AccessCache;
  cache->Entry = entry;
  cache->Offset = entry->Multiplicity - 1;
  
  for ( unsigned int idx = 1; idx < this->DataSize; ++idx ) 
    {
    if ( (entry->Value == this->Data[idx]) &&  // still same data value
	 (entry->Multiplicity < 255) ) 
      { // and no counter overrun
      ++entry->Multiplicity;
      } 
    else 
      {
      ++entry; // otherwise: start new RLE entry
      entry->Value = this->Data[idx];
      entry->Multiplicity = 1;
      }
    // original index multiple of CacheBlockLength?
    // if so, create next access cache entry.
    if ( !(idx % CacheBlockLength) ) 
      {
      // yes: create next cache table entry
      ++cache;
      cache->Entry = entry;
      cache->Offset = entry->Multiplicity - 1;
      }
    }
  // add guard entry (item count == zero) at the end of encoded sequence
  ++entry;
  entry->Multiplicity = 0;
  ++cache;
  cache->Entry = NULL;
  
  // did we have to convert original array data?
  if ( dataArray->GetType() != this->GetType() ) 
    {
    // yes: delete intermediate data array.
    this->Free( this->Data );
    }
  this->Data = NULL;
}

template<class T>
void TemplateArrayRLE<T>::GetRaw
( T *const values, const int index, const unsigned length ) const 
{
  CheckBounds( index+length-1, this->DataSize );
  const typename Self::AccessCacheEntry *cache = AccessCache + (index / CacheBlockLength);
  const typename Self::RLEentry *entry = cache->Entry;

  // number of data values to skip from beginning of cache entry block
  unsigned int skip = index % CacheBlockLength + cache->Offset;
  
  // skip as many RLE encoded bytes as necessary
  while ( skip >= entry->Multiplicity ) 
    {
    skip -= entry->Multiplicity;
    ++entry;
    }
  
  T *nextValue = values;
  unsigned toGo = length;
  unsigned remain = entry->Multiplicity - skip;
  while ( toGo ) 
    {
    unsigned fromThisEntry = std::min( toGo, remain );
    
    if ( this->PaddingFlag && ( this->Padding == entry->Value)) 
      {
      memset( nextValue, 0, fromThisEntry * sizeof( *nextValue ) );
      nextValue += fromThisEntry;
      } 
    else
      {
      for ( unsigned idx = 0 ; idx < fromThisEntry; ++idx, ++nextValue ) 
	{
	*nextValue = entry->Value;
	}
      }
    
    toGo -= fromThisEntry;
    ++entry;
    remain = entry->Multiplicity;
    }
}

template<class T>
bool TemplateArrayRLE<T>::GetRangeTemplate ( T& min, T& max ) const 
{
  // find first non-Padding element
  const typename Self::RLEentry *entry = DataRLE;
  if ( this->PaddingFlag ) 
    while ( entry->Multiplicity && (entry->Value == this->Padding) )
      ++entry;
  
  // didn't find any? return with error flag.
  if ( ! entry->Multiplicity ) 
    {
    return false;
    }
  
  // otherwise: search for min/max from here
  min = entry->Value;
  max = min;

  if ( this->PaddingFlag ) 
    {
    for ( ; entry->Multiplicity; ++entry ) 
      {
      if ( entry->Value != this->Padding ) 
	{
	if (entry->Value > max) max = entry->Value;
	if (entry->Value < min) min = entry->Value;
	}
      }
    } 
  else
    {
    for ( ; entry->Multiplicity; ++entry ) 
      {
      if (entry->Value > max) max = entry->Value;
      if (entry->Value < min) min = entry->Value;
      }
    }
  return true;
}

template<class T>
Histogram<unsigned int>*
TemplateArrayRLE<T>::GetHistogram( const unsigned int numberOfBins ) const
{
  // find first non-Padding element
  const typename Self::RLEentry *entry = DataRLE;
  if ( this->PaddingFlag ) 
    while ( entry->Multiplicity && (entry->Value == this->Padding) )
    ++entry;

  // didn't find any? return with error value (NULL).
  if ( ! entry->Multiplicity ) 
    {
    return NULL;
    }
  
  // otherwise: create and fill histogram
  Histogram<unsigned int>* histogram = new Histogram<unsigned int>( numberOfBins );
  T min, max;
  this->GetRangeTemplate( min, max );
  histogram->SetRange( min, max );
  
  min = entry->Value;
  max = min;
  
  if ( this->PaddingFlag ) 
    {
    for ( ; entry->Multiplicity; ++entry ) 
      {
      if ( entry->Value != this->Padding ) 
	{
	histogram->Increment( histogram->ValueToBin( entry->Value ), entry->Multiplicity );
	}
      }
    } 
  else
    {
    histogram->Increment( histogram->ValueToBin( entry->Value ), entry->Multiplicity );
    }
  
  return histogram;
}

template<class T>
unsigned TemplateArrayRLE<T>::RequiredLengthRLE() const
{
  unsigned requiredLength = 1;
  T currentData = *this->Data;
  
  byte multiplicity = 1;
  for ( unsigned int idx = 1; idx < this->DataSize; ++idx, ++multiplicity ) 
    {
    if ( (this->Data[idx] != currentData) || (multiplicity==255) ) 
      {
      ++requiredLength;
      currentData = this->Data[idx];
      multiplicity = 0;
      }
    }
  
  return requiredLength;
}

/// Array of (unsigned) byte values with Run Length Encoding.
template class TemplateArrayRLE<byte>;

/// Array of (signed) char values with Run Length Encoding.
template class TemplateArrayRLE<char>;

/// Array of signed short values with Run Length Encoding.
template class TemplateArrayRLE<short>;

/// Array of unsigned short values with Run Length Encoding.
template class TemplateArrayRLE<unsigned short>;

/// Array of (signed) integer values with Run Length Encoding.
template class TemplateArrayRLE<int>;

/// Array of single-precision float values with Run Length Encoding.
template class TemplateArrayRLE<float>;

/// Array of double-precision float values with Run Length Encoding.
template class TemplateArrayRLE<double>;

} // namespace cmtk
