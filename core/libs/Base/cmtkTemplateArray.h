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

#ifndef __cmtkTemplateArray_h_included_
#define __cmtkTemplateArray_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypedArray.h>
#include <Base/cmtkTypedArrayFunction.h>

#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Template for Variable-Typed Data Arrays.
 * From this object, various children are derived for the concrete data types
 * to be handled.
 * @author Torsten Rohlfing 
 */
template<class T>
class TemplateArray : 
  /** Inherit class interface from virtual base class. */
  public TypedArray
{
public:
  /// This type.
  typedef TemplateArray<T> Self;

  /// Base class.
  typedef TypedArray Superclass;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Data type traits. */
  typedef DataTypeTraits<T> TypeTraits;

  /// Create array of this type.
  static typename Self::SmartPtr Create( const size_t size )
  {
    return typename Self::SmartPtr( new Self( size ) );
  }

  /// Return const pointer to actual data.
  const T* GetDataPtrConcrete() const { return this->Data; }

  /// Return pointer to actual data.
  T* GetDataPtrConcrete() { return this->Data; }

  /// Set the value to mark non-existent data using template data type.
  void SetPaddingValueTemplate ( const T paddingData ) 
  {
    this->Padding = paddingData;
    this->PaddingFlag = true;
  }
  
  /** Constructor.
   * A typed array is built from an existing array.
   *@param data Pointer to the array of values to be stored.
   *@param datasize Number of elements in the data array.
   *@param freeArray Tag of the memory handler that has to be used for freeing
   * the array after use.
   *@param paddingflag Flag that indicates whether there are missing elements in
   * the existing data array.
   *\param paddingData Value that marks missing elements in the data array if "paddingFlag" is true.
   */
  TemplateArray ( void *const data, const size_t datasize, const bool freeArray, const bool paddingflag, const void* paddingData ) 
  {
    m_DataType = TypeTraits::GetScalarDataType();
    Data = static_cast<T*>( data );
    DataSize = datasize;
    FreeArray = freeArray;
    PaddingFlag = paddingflag;
    if ( paddingData )
      Padding = *((T*)paddingData);
    else
      Padding = (T) 0;
  }

  /** Constructor.
   * Allocate an array of a given size.
   */
  TemplateArray ( const size_t datasize = 0 ) 
  {
    m_DataType = TypeTraits::GetScalarDataType();
    Data = NULL;
    this->Alloc( datasize ); 
  }

  /** Destructor.
   * Free memory by a call to FreeData().
   */
  virtual ~TemplateArray () 
  {    
    this->FreeData();
  }
  
  /// Set the value to mark non-existent data using interface data type.
  virtual void SetPaddingValue ( const Types::DataItem paddingData ) 
  {
    this->SetPaddingValueTemplate( TypeTraits::Convert( paddingData ) );
  }

  /// Set the value to mark non-existent data using interface data type.
  virtual void SetPaddingPtr ( const void* paddingData ) 
  {
    this->SetPaddingValueTemplate( *((T*) paddingData) );
  }

  /// Replace Padding data (padding) with given value.
  virtual void ReplacePaddingData ( const Types::DataItem value = 0 ) 
  {
    if ( PaddingFlag ) 
      {
      T v = TypeTraits::Convert( value );
      for ( size_t i = 0; i < DataSize; ++i )
	if ( Data[i] == Padding ) 
	  {
	  Data[i] = v;
	  }
      }
  }
  
  /** Check for NULL data at given index.
   * If this array does not have PaddingFlag set, the result is always false.
   */
  virtual bool IsPaddingAt( const size_t index ) const 
  {
    return PaddingFlag && (Data[index]==Padding);
  }
  
  /** Check for PADDING data or zero value at given index.
   */
  virtual bool IsPaddingOrZeroAt( const size_t index ) const 
  {
    return (PaddingFlag && (Data[index]==Padding)) || Data[index] == (T)0;
  }
  
  virtual void SetPaddingAt ( const size_t index = 0 ) 
  {
    if ( !PaddingFlag ) 
      {
      Padding = TypeTraits::ChoosePaddingValue();
      PaddingFlag = true;
      }
    Data[index] = Padding;
  }
  
  /// Return id of primitive data type handled by this object.
  virtual ScalarDataType GetType () const { return TypeTraits::GetScalarDataType(); }

  /// Return size in bytes of the primitive data type handled by this object.
  virtual size_t GetItemSize () const { return sizeof(T); }

  /// Return pointer to an element in this objects data array.
  virtual void* GetDataPtr( const size_t offset = 0 ) 
  {
    return Data + offset; 
  }

  /// Return const pointer to an element in this objects data array.
  virtual const void* GetDataPtr( const size_t offset = 0 ) const 
  {
    return Data + offset; 
  }
  
  /// Return pointer to an element in this objects data array.
  virtual T* GetDataPtrTemplate( const size_t offset = 0 ) 
  {
    return Data + offset; 
  }
  
  /// Return const pointer to an element in this objects data array.
  virtual const T* GetDataPtrTemplate( const size_t offset = 0 ) const 
  {
    return Data + offset; 
  }
  
  /** Convert and copy continuous sub-array to a given destination.
   *@param toPtr A pointer to the location where the data shall be stored.
   *@param fromIdx The index of the first copied data item in the array,
   * beginning with 0.
   *@param len Length, ie. number of values, to copy. The calling function
   * must take care that there is enough memory pointed to by toPtr to store
   * this many values of the transfer data type.
   *@param substPadding Where there is padding data in the copied range, this value
   * is put to the output.
   *@return The pointer given to this function to store the desired data to.
   */
  virtual Types::DataItem* GetSubArray( Types::DataItem *const toPtr, const size_t fromIdx , const size_t len, const Types::DataItem substPadding = 0 ) const 
  {
    int index = fromIdx;
    
    if ( PaddingFlag ) 
      {
      for ( size_t i = 0; i<len; ++i, ++index ) 
	{
	T value = Data[index];
	if (value == Padding)
	  toPtr[i] = substPadding;
	else
	  toPtr[i] = static_cast<Types::DataItem>( value );
	}
      } 
    else
      {
      for ( size_t i = 0; i<len; ++i, ++index )
	toPtr[i] = static_cast<Types::DataItem>( Data[index] );
      }
    return toPtr;
  }

  /** Return a sub array of the current instance.
   * The data is returned in a newly allocated primitive array of the Types::DataItem
   * data exchange type.
   *@param fromIdx Copying starts at this index in the array. The range for
   * this parameter is [0..Size-1].
   *@param len Number of data elements to be copied.
   *@param substPadding If this flag is set (default is NO), then during copying
   * all elements marked as "non-existent" are replaced by zero.
   *@return Pointer to a newly allocated Types::DataItem array. Allocation is done
   * using the allocator given as template parameter "M" to this class.
   */
  virtual Types::DataItem* GetSubArray( const size_t fromIdx, const size_t len, const Types::DataItem substPadding = 0 ) const 
  {
    Types::DataItem* buffer = Memory::AllocateArray<Types::DataItem>( len );
    return this->GetSubArray( buffer, fromIdx, len, substPadding );
  }

  /** Convert Types::DataItem to template type.
   * This function takes an Types::DataItem value and converts it into a value of the
   * type stored in this array.
   *@param value The item in Types::DataItem representation.
   *@return The value of type T corresponding to the "value" parameter's value.
   * Conversion is done using the type converter class given as template
   * parameter "C" to this class. Therefore, rounding will occur if necessary.
   */
  virtual T ConvertItem ( const Types::DataItem value ) 
  {
    return TypeTraits::Convert( value, PaddingFlag, Padding );
  }

  /** Convert to typed array of any given template type.
   */
  virtual TypedArray::SmartPtr Convert(  const ScalarDataType dtype ) const;


  /** Convert a sub-array to any given primitive data type.
   *\todo It would probably be a good idea to preserve PaddingData as much as
   * possible during the conversion.
   */
  virtual void* ConvertSubArray( const ScalarDataType dtype, const size_t fromIdx, const size_t len ) const;
  
  /** Convert a sub-array to given primitive data type into existing array.
   */
  virtual void ConvertSubArray( void *const destination, const ScalarDataType dtype, const size_t fromIdx, const size_t len ) 
    const;

  /** Change endianness of data.
   */
  virtual void ChangeEndianness();

  /** Scale values in the array.
   * A call to this member function will perform an in-place rescaling of the
   * values in the array.
   *@param scale The original data value is multiplied by the parameter first.
   *@param offset This value is added to the original data value after 
   * multiplying it by the scale parameter.
   */
  virtual void Rescale( const Types::DataItem scale = 1, const Types::DataItem offset = 0 ) 
  {
#pragma omp parallel for    
    for ( size_t i = 0; i < this->DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) )
	Data[i] = TypeTraits::Convert( (scale * Data[i]) + offset );
  }
  
  /** Scale and shift values in the array.
   * Data values are scaled first, then the offsewt is added, and finally the (left) bit shift is applied.
   * In fact, to make sure we're not messing up floats, the bit shift is applied as a multiplication.
   */
  virtual void RescaleAndShift( const Types::DataItem scale = 1, const Types::DataItem offset = 0, const size_t shiftBits = 0 ) 
  {
    const long int shiftMultiplier = (1<<shiftBits);
#pragma omp parallel for    
    for ( size_t i = 0; i < this->DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) )
	Data[i] = TypeTraits::Convert( ((scale * Data[i]) + offset) * shiftMultiplier );
  }
  
  /** Scale values in the array with truncation boundaries.
   * A call to this member function will perform an in-place rescaling of the
   * values in the array with value range truncation. Truncation takes place
   * after the scaling itself, i.e., the truncation boundaries refer to the
   * already scaled values.
   *@param scale The original data value is multiplied by the parameter first.
   *@param offset This value is added to the original data value after 
   * multiplying it by the scale parameter.
   *@param truncLo Lower truncation boundary. Scaled items below this 
   * threshold will be set to equal its value.
   *@param truncHi Upper truncation boundary. Scaled items above this 
   * threshold will be set to equal its value.
   */
  virtual void Rescale( const Types::DataItem scale, const Types::DataItem offset, const Types::DataItem truncLo, const Types::DataItem truncHi = CMTK_ITEM_MAX ) 
  {
#pragma omp parallel for    
    for ( size_t i = 0; i < this->DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) ) 
	{
	Data[i] = TypeTraits::Convert( (scale * Data[i]) + offset );
	if ( Data[i] < truncLo )
	  Data[i] = TypeTraits::Convert( truncLo );
	else
	  if ( Data[i] > truncHi )
	    Data[i] = TypeTraits::Convert( truncHi );
	}
  }

  /** Apply gamma correction.
   *@param gamma The gamma correction coefficient.
   */
  virtual void GammaCorrection( const Types::DataItem gamma );

  /** Apply real function to data.
   */
  virtual void ApplyFunctionFloat( typename Self::FunctionTypeFloat f );

  /** Apply real function to data.
   */
  virtual void ApplyFunctionDouble( typename Self::FunctionTypeDouble f );

  /** Threshold data.
   * All values above upper threshold are set to upper thrershold. All values
   * below lower threshold are set to lower threshold.
   */
  virtual void Threshold( const Types::DataItemRange& range ) 
  {
    const T lo = TypeTraits::Convert( range.m_LowerBound );
    const T hi = TypeTraits::Convert( range.m_UpperBound );
#pragma omp parallel for    
    for ( size_t i = 0; i < this->DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) ) 
	{
	if ( Data[i] < lo )
	  Data[i] = lo;
	else
	  if ( Data[i] > hi )
	    Data[i] = hi;
	}
  }
  
  /** Threshold data.
   * All values outside the threshold range are set to the Padding (padding)
   * value.
   */
  virtual void ThresholdToPadding( const Types::DataItemRange& range ) 
  {
    const T lo = TypeTraits::Convert( range.m_LowerBound );
    const T hi = TypeTraits::Convert( range.m_UpperBound );
#pragma omp parallel for    
    for ( size_t i = 0; i < this->DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) ) 
	{
	if ( (Data[i] < lo) || (Data[i] > hi) )
	  Data[i] = this->Padding;
	}
  }
  
  /** Binarize array values with given threshold.
   * All values greater than threshold (default: zero) are set to one, all
   * values smaller or equal are set to zero.
   */
  virtual void Binarize( const Types::DataItem threshold = 0 ) 
  {
    T thresh = TypeTraits::Convert( threshold );
    
    T one = TypeTraits::Convert( 1.0 ), zero = TypeTraits::Convert( 0.0 );
#pragma omp parallel for    
    for ( size_t i = 0; i < this->DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) ) 
	{
	if ( Data[i] > thresh )
	  Data[i] = one;
	else
	  Data[i] = zero;
	}
  }
  
  /// Convert all values to absolute values.
  virtual void MakeAbsolute()
  {
#pragma omp parallel for    
    for ( size_t i = 0; i < this->DataSize; ++i )
      if ( ! PaddingFlag || (Data[i] != Padding ) )
	Data[i] = TypeTraits::Abs( Data[i] );
  }
  
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
    CheckBounds( index, DataSize );
    if (PaddingFlag && (Padding==Data[index])) 
      {
      value = 0;
      return false;
      }
    value = static_cast<Types::DataItem>( Data[index] );
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
    CheckBounds( index+length-1, DataSize );
    for ( size_t i = 0; i < index+length; ++i )
      if (PaddingFlag && (Padding==Data[index]))
	values[i] = 0;
      else
	values[i] = static_cast<Types::DataItem>( Data[index] );
  }

  /** Set an item in the array.
   *@param value The new value for the specified item.
   *@param index The index of the item to set. Valid values are in the 
   * range [0..Datasize()-1].
   *@see Convert
   */
  virtual void Set ( const Types::DataItem value, const size_t index ) 
  {
    CheckBounds( index, DataSize );
    Data[index] = this->ConvertItem( value );
  }

  /** Return padding data value as pointer to native representation.
   * This function returns a pointer to the value used by this object for 
   * marking non-existent data. As the type of this value depends on the
   * template parameters of this class, only a void pointer can be returned.
   * The caller has to interpret the value stored at that location itself.
   *
   * If this array does NOT have a padding data value, the data pointed to by the
   * result of this function is undefined.
   *@return Pointer to the padding value of this object.
   */
  virtual void* GetPaddingPtr () const { return (void*)&Padding; }

  /** Return padding data value as pointer to native representation.
   */
  virtual Types::DataItem GetPaddingValue() const { return static_cast<Types::DataItem>( Padding ); }

  /// Test if there is PADDING data at a particular location.
  virtual bool PaddingDataAt ( const size_t index ) const 
  {
    return PaddingFlag && (Data[index] == Padding);
  }
  
  /** Get the whole array data as an exchange type array.
   *@return Pointer to a memory region allocated by Memory::AllocateArray(). This region is
   * filled with all values in the present array as Types::DataItem values. The created
   * array is not maintained by this object. The caller has to make sure free()
   * is called for it.
   */
  virtual Types::DataItem* GetData () const 
  {
    Types::DataItem* Result = Memory::AllocateArray<Types::DataItem>( DataSize );
    if ( Result ) 
      {
      for ( size_t idx = 0; idx < DataSize; ++idx )
	Result[idx] = (Types::DataItem) Data[idx];
      }
    return Result;
  }

  /** Set all data from an Types::DataItem array.
   * This function sets all values stored in the present array from a memory
   * region with Types::DataItem values.
   *@param data Pointer to an array of Types::DataItem values.
   * Control over the source array is not taken by this object. If it is on the heap, 
   * then the calling routine remains responsible for de-allocating the array afterwards.
   */
  virtual void SetData ( Types::DataItem *const data ) 
  {
#pragma omp parallel for    
    for ( size_t idx = 0; idx < this->DataSize; ++idx )
      Data[idx] = this->ConvertItem( data[idx] );
  }

  /** Clear entire array.
   *@param usePaddingData If this flag is set, then the array will be filled with
   * the PaddingData value, if one exists. Otherwise, the array will be filled
   * with the respective data type's zero value.
   */
  virtual void ClearArray ( const bool usePaddingData = false ) 
  {
    if ( usePaddingData && PaddingFlag ) 
      {
      for ( size_t idx = 0; idx < DataSize; ++idx )
	Data[idx] = Padding;
      } 
    else
      {
      memset( Data, 0, sizeof( *Data ) * this->GetDataSize() );
      }
  }

  /// Calculate minimum and maximum data value.
  virtual const Types::DataItemRange GetRange() const;

  /// Calculate minimum and maximum data value.
  virtual const Types::Range<T> GetRangeTemplate() const;

  /// Calculate entropy of distribution of values in this array.
  virtual double GetEntropy( const bool fractional = CMTK_HISTOGRAM_DISCRETE, const int numberOfBins = 128 ) const;
  
  virtual double GetEntropy( Histogram<unsigned int>& histogram ) const;
  virtual double GetEntropy( Histogram<double>& histogram, const bool fractional = CMTK_HISTOGRAM_DISCRETE ) const;
  virtual double GetEntropy( Histogram<double>& histogram, const double* kernel, const size_t kernelRadius ) const;

  /** Calculate statistics.
   * Results will be both zero if there is not data in the array.
   *@return The number of valid (i.e., non-padding) values that constitute the
   * given results.
   */
  virtual size_t GetStatistics ( Types::DataItem& mean, Types::DataItem& variance ) const;

  /** Set data block to constant value.
   */
  virtual void BlockSet( const Types::DataItem value, const size_t fromOffset, const size_t toOffset );

  /** Get data histogram.
   *@return A histogram object filled with the relative frequencies of values 
   * in this array.
   */
  virtual Histogram<unsigned int>::SmartPtr GetHistogram( const unsigned int numberOfBins /*!< Number of histogram bins */,
							  const bool centeredBins = false /*!< Flag for bins centered around the samples*/ ) const;

  virtual void ApplyFunctionObject( const TypedArrayFunction& f );

protected:
  /// Clone this object.
  virtual Self* CloneVirtual() const
  {
    Self* clone = new Self( this->DataSize );
    memcpy( clone->Data, this->Data, this->DataSize * sizeof( T ) );
    clone->Padding = this->Padding;
    clone->PaddingFlag = this->PaddingFlag;
    clone->m_DataClass = this->m_DataClass;
    return clone;
  }

private:
  /// The acutal data array.
  T *Data;

  /// Value used for missing data.
  T Padding;

  /** Allocate data array.
   *@param datasize Number of data items to allocate memory for.
   */
  virtual void Alloc ( const size_t datasize ) 
  {
    DataSize = datasize;
    if ( DataSize ) 
      {
      if ( Data ) 
	{
	Memory::DeleteArray( Data );
	}
      Data = Memory::AllocateArray<T>( DataSize );
      if ( Data == NULL ) 
	{
	this->DataSize = 0;
	}
      FreeArray = true;
      } 
    else
      {
      Data = NULL;
      FreeArray = false;
      }
  }
  
  /** De-allocate data array.
   * The array is freed using the same memory handler that allocated it.
   * In any case, the current Data pointer is set to NULL.
   */
  virtual void FreeData() 
  {
    if ( Data && FreeArray ) 
      {
      Memory::DeleteArray( Data );
      } 
    Data = NULL;
  }
};

/**@name Shortcut class typedefs for typed arrays. */
//@{

/// Array of (unsigned) byte values.
typedef TemplateArray<byte>   ByteArray;

/// Array of (signed) char values.
typedef TemplateArray<char>   CharArray;

/// Array of signed short values.
typedef TemplateArray<short>  ShortArray;

/// Array of unsigned short values.
typedef TemplateArray<unsigned short> UShortArray;

/// Array of (signed) integer values.
typedef TemplateArray<int>    IntArray;

/// Array of single-precision float values.
typedef TemplateArray<float>  FloatArray;

/// Array of double-precision float values.
typedef TemplateArray<double> DoubleArray;

//@}

//@}

} // namespace cmtk

#include "cmtkTemplateArray.txx"
#include "cmtkTemplateArray_Statistics.txx"

#endif // #ifndef __cmtkTemplateArray_h_included_
