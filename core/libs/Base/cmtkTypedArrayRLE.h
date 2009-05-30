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

#ifndef __cmtkTypedArrayRLE_h_included_
#define __cmtkTypedArrayRLE_h_included_

#include <cmtkconfig.h>

#include <cmtkTypedArray.h>
#include <cmtkTemplateArray.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Replacement for TypedArray classes using Run Length Encoding.
 * This class is intended to handle scalar data with large blocks of repeating
 * values. It is, for example, particularly well suited for storing large
 * label fields.
 */
template<class T>
class TemplateArrayRLE :
  /// Inherit some stuff from usual typed array.
  public TemplateArray<T>
{
public:
  /** This class. */
  typedef TemplateArrayRLE<T> Self;

  /** Smart pointer to this class. */
  typedef SmartPointer<Self> SmartPtr;

  /** Superclass. */
  typedef TemplateArray<T> Superclass;

  /** Data type traits. */
  typedef typename Superclass::TypeTraits TypeTraits;

  /** Constructor.
   *@param dataArray Array of values that will be RLE encoded and make up the
   * newly constructed array.
   *@param cacheBlockLength The block length for the internal rapid-access
   * cache. This can be provided by a client to match its requirements. For
   * example, when encoding a 2-D or 3-D image that is to be read only by
   * complete rows, this could be the number of samples per row. If this
   * parameter is not given, the block length is taken from the constant
   * DefaultCacheBlockLength.
   */
  TemplateArrayRLE( const TypedArray* dataArray, const unsigned cacheBlockLength = 0 );

  /// Destructor.
  virtual ~TemplateArrayRLE() 
  {
    delete[] AccessCache;
    delete[] DataRLE;
  }

  /** Return the cache block length used by this object.
   * Clients of this class can make reading of data in blocks more efficient
   * if they read chunks of the given length the offsets of which are multiples
   * of the block length.
   */
  unsigned GetCacheBlockLength() const { return CacheBlockLength; }

  /** Get an item from the array.
   *@param value The variable referenced by this parameter is set to the
   * data item stored at the given location in the array. If this item is
   * marked is "padding data", ie. non-existent, "value" is set to zero.
   *@param index The index of the item to retrieve. Valid values are in the 
   * range [0..Datasize()-1].
   *@return A non-zero value is returned if and only if the value stored in
   * the array at the given location is marked as valid data.
   */
  virtual bool Get ( Types::DataItem& value, const size_t index ) const 
  {
    CheckBounds( index, this->DataSize );
    const typename Self::AccessCacheEntry *cache = AccessCache + (index / Self::CacheBlockLength);
    // number of data values to skip from beginning of cache entry block
    unsigned int remain = index % Self::CacheBlockLength + cache->Offset;
    
    const typename Self::RLEentry *entry = cache->Entry;
    // skip as many RLE encoded bytes as necessary
    while ( remain > entry->Multiplicity ) 
      {
      remain -= entry->Multiplicity;
      ++entry;
      }
    if ( this->PaddingFlag && ( this->Padding == entry->Value)) 
      {
      value = 0;
      return false;
      }
    value=(Types::DataItem) entry->Value;
    return true;
  }

  /** Get a sequence of items from the array.
   *@param values This must point to an allocated array of at least as many
   * Types::DataItem objects as given in the "length" parameter.
   *@param index The index of the item to retrieve. Valid values are in the 
   * range [0..Datasize()-1].
   *@param length Number of consecutive values to get.
   */
  virtual void GetSequence ( Types::DataItem *const values, const size_t index, const size_t length ) const
  {
    CheckBounds( index+length-1, this->DataSize );
    const typename Self::AccessCacheEntry *cache = AccessCache + (index / Self::CacheBlockLength);
    // number of data values to skip from beginning of cache entry block
    unsigned int remain = index % Self::CacheBlockLength + cache->Offset;
    
    const typename Self::RLEentry *entry = cache->Entry;
    // skip as many RLE encoded bytes as necessary
    while ( remain > entry->Multiplicity ) 
      {
      remain -= entry->Multiplicity;
      ++entry;
      }
    
    Types::DataItem *nextValue = values;
    unsigned toGo = length;
    while ( toGo ) 
      {
      if ( this->PaddingFlag && ( this->Padding == entry->Value)) 
	{
	*nextValue = 0;
	}
      *nextValue=(Types::DataItem) entry->Value;
      ++nextValue;
      --toGo;
      
      // decrement number of values left in current RLE entry
      --remain;
      // is entry finished?
      if ( ! remain ) 
	{
	// move to next RLE field
	++entry;
	remain = entry->Multiplicity;
	}
      }
  }
  
  /** Get a sequence of items from the array in this objects native type.
   *@param values This must point to an allocated array of at least as many
   * Types::DataItem objects as given in the "length" parameter.
   *@param index The index of the item to retrieve. Valid values are in the 
   * range [0..Datasize()-1].
   *@param length Number of consecutive values to get.
   */
  virtual void GetRaw ( T *const values, const int index, const unsigned length ) const;
  
  /** Calculate minimum and maximum data value.
   * We can do this very efficiently for RLE data, since each RLE block only
   * needs to be considered once, rather than on a value-by-value basis.
   */
  virtual bool GetRangeTemplate( T& min, T& max ) const;
  
  /// Convert all values to absolute values.
  virtual void MakeAbsolute() 
  {
    if ( this->PaddingFlag ) 
      {
      for ( typename Self::RLEentry *entry = DataRLE; entry->Multiplicity; ++entry ) 
	{
	if ( entry->Value != this->Padding )
	  entry->Value = TypeTraits::Abs( entry->Value );
	}
      } 
    else
      {
      for ( typename Self::RLEentry *entry = DataRLE; entry->Multiplicity; ++entry ) 
	{
	entry->Value = TypeTraits::Abs( entry->Value );
	}
      }
  }
  
  /** Get data histogram.
   *@return A histogram object filled with the relative frequencies of values 
   * in this array.
   */
  virtual Histogram<unsigned int>* GetHistogram (const unsigned int numberOfBins ) const;

private:
  /** Default length of blocks in the access cache.
   * Smaller numbers will allow for faster access, but result in a larger
   * portion of memory required to store the cache table.
   */
  static const unsigned DefaultCacheBlockLength = 256;

  /// Entry in the RLE encoded array
  class RLEentry 
  {
  public:
    /// Number of successive occurences of this value.
    byte Multiplicity;
    /// Actual data value.
    T Value;
  };
  
  /// RLE encoded data.
  RLEentry *DataRLE;

  /// Number of RLE entries.
  unsigned NumberEntriesRLE;

  /// Entry for the fast access cache table.
  class AccessCacheEntry 
  {
  public:
    /// Index of the RLE entry containing the given data item.
    typename Self::RLEentry *Entry;

    /// Offset inside the block described by the RLE entry.
    byte Offset;
  };

  /// Cache table for fast access to explicitly indexed data items.
  AccessCacheEntry *AccessCache;

  /// Block length in the access cache.
  unsigned CacheBlockLength;

  /// Compute number of RLE entries necessary to store given uncompressed data.
  unsigned RequiredLengthRLE() const;

};

/// Array of (unsigned) byte values with Run Length Encoding.
typedef TemplateArrayRLE<byte> ByteArrayRLE;

/// Array of (signed) char values with Run Length Encoding.
typedef TemplateArrayRLE<char> CharArrayRLE;

/// Array of signed short values with Run Length Encoding.
typedef TemplateArrayRLE<short> ShortArrayRLE;

/// Array of unsigned short values with Run Length Encoding.
typedef TemplateArrayRLE<unsigned short> UShortArrayRLE;

/// Array of (signed) integer values with Run Length Encoding.
typedef TemplateArrayRLE<int> IntArrayRLE;

/// Array of single-precision float values with Run Length Encoding.
typedef TemplateArrayRLE<float> FloatArrayRLE;

/// Array of double-precision float values with Run Length Encoding.
typedef TemplateArrayRLE<double> DoubleArrayRLE;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypedArrayRLE_h_included_
