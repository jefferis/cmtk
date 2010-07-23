/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#ifndef __cmtkTypedArray_h_included_
#define __cmtkTypedArray_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkMacros.h"
#include "Base/cmtkTypes.h"
#include "Base/cmtkMathUtil.h"
#include "Base/cmtkDataTypeTraits.h"
#include "Base/cmtkHistogram.h"
#include "Base/cmtkTypedArrayFunction.h"

#include "System/cmtkSmartPtr.h"
#include "System/cmtkSmartConstPtr.h"
#include "System/cmtkException.h"

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <climits>

#include <memory.h>

#include <vector>

#ifdef DEBUG
#define CheckBounds(index,bound) \
  if (!(index<bound)) throw( Exception( "Index is outside bounds", this ) );
#else
#define CheckBounds(index,bound)
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Generic Variable-Typed Data Array.
 * This class delivers the common interface of all variable-typed data arrays.
 *@author Torsten Rohlfing 
 */
class TypedArray
{
  /// Parameter defining the primitive data class.
  cmtkGetSetMacro(DataClass,DataClass);

public:
  /// This class.
  typedef TypedArray Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /** Create typed data array from existing array of values.
   *@param dtype Type specifier.
   *@param data Pointer to the existing data array.
   *@param size Number of data items in the array.
   *@param freeArray On destruction, the created object will de-allocate the
   * array's memory if and only if this flag is true.
   *@param paddingflag If this flag is not zero, padding data exists in the array.
   *@param padding Value used for padding data.
   *@return A pointer to a new typed array object, or NULL if an error 
   * occurred.
   */
  static Self::SmartPtr Create
  ( const ScalarDataType dtype, void *const data, const size_t size, const bool freeArray = true, const bool paddingFlag = false, const void* paddingData = NULL );
  
  /** Create typed data array, allocating new memory for the items array.
   *@param dtype Type specifier.
   *@param size Number of items in the array to be allocated. Memory will be
   * freed when the created object is destroyed. The values in the array are 
   * not initialized.
   *@return A pointer to a new typed array object, or NULL if an error
   * occurred.
   */
  static Self::SmartPtr Create( const ScalarDataType dtype, const size_t size );

  /** Free data allocated by this object. */
  static void Free( void *const data )
  {
    free( data );
  }
  
  /** Get an item from the specified index in the data array.
   * If there is no valid data present at the given index, zero is returned
   * and the resulting item is zero, too. This can only happen if PaddingFlag
   * is not zero.
   */
  virtual bool Get ( Types::DataItem&, const size_t ) const = 0;

  /** Get a sequence of items from the array.
   *@param values This must point to an allocated array of at least as many
   * Types::DataItem objects as given in the "length" parameter.
   *@param index The index of the item to retrieve. Valid values are in the 
   * range [0..Datasize()-1].
   *@param length Number of consecutive values to get.
   */
  virtual void GetSequence ( Types::DataItem *const values, const size_t index, const size_t length ) const = 0;

  /// Sets the specified array element.
  virtual void Set ( const Types::DataItem, const size_t ) = 0;

  /// Return type of stored data.
  virtual ScalarDataType GetType () const = 0;

  /// Return the number of bytes per stored item.
  virtual size_t GetItemSize () const = 0;

  /** Get the adress of the real value used for no data present.
   */
  virtual void* GetPaddingPtr () const = 0;

  /// Get the transfer representation of the value used for no data present.
  virtual Types::DataItem GetPaddingValue () const = 0;

  /// Test if there is NULL data at a particular location.
  virtual bool PaddingDataAt ( const size_t index ) const = 0;

  /** Get complete aray of data values.
   * The calling routine is responsible for de-allocating the returned array
   * by calling 'free' after use.
   */
  virtual Types::DataItem* GetData () const = 0;

  /** Set all data from an Types::DataItem array.
   * This function sets all values stored in the present array from a memory
   * region with Types::DataItem values.
   *@param data Pointer to an array of Types::DataItem values.
   * Control over the source array is not taken by this object. If it is on the heap, 
   * then the calling routine remains responsible for de-allocating the array afterwards.
   */
  virtual void SetData( Types::DataItem *const data ) = 0;

  /** Convert to typed array of any given template type.
   */
  virtual Self::SmartPtr Convert(  const ScalarDataType dtype ) const = 0;

  /** Convert a sub-array to any given primitive data type.
   *\return Newly allocated memory of given type. Caller is responsible for freeing this memory with a call to
   * Memory::FreeArray().
   */
  virtual void* ConvertSubArray( const ScalarDataType dtype, const size_t fromIdx, const size_t len ) const = 0;

  /** Convert a sub-array to given primitive data type into existing array.
   */
  virtual void ConvertSubArray( void *const destination, const ScalarDataType dtype, const size_t fromIdx, const size_t len ) const = 0;

  /** Convert the array to any given data type.
   * This function uses ConvertSubArray to convert the complete array.
   *\return Newly allocated memory of given type. Caller is responsible for freeing this memory with a call to
   * Memory::FreeArray().
   *@see ConvertSubArray
   */
  virtual void* ConvertArray ( const ScalarDataType dtype ) const 
  {
    return this->ConvertSubArray( dtype, 0, DataSize );
  }
  
  /** Clear entire array.
   * This method is implemented by derived template classes for better
   * efficiency.
   *@param usePaddingData If this flag is set, then the array will be filled with
   * the PaddingData value, if one exists. Otherwise, the array will be filled
   * with the respective data type's zero value.
   */
  virtual void ClearArray ( const bool usePaddingData = false ) = 0;

  /** Change endianness of data.
   */
  virtual void ChangeEndianness() = 0;

  /** Scale values in the array.
   * A call to this member function will perform an in-place rescaling of the
   * values in the array.
   *@param scale The original data value is multiplied by the parameter first.
   *@param offset This value is added to the original data value after 
   * multiplying it by the scale parameter.
   */
  virtual void Rescale( const Types::DataItem scale = 1, const Types::DataItem offset = 0 ) = 0;

  /** Scale and shift values in the array.
   */
  virtual void RescaleAndShift( const Types::DataItem scale = 1, const Types::DataItem offset = 0, const size_t shiftBits = 0 ) = 0;

  /** Scale values in the array to a given target range.
   */
  virtual void RescaleToRange( const Types::DataItemRange& toRange );

  /** Apply gamma correction.
   *@param gamma The gamma correction coefficient.
   */
  virtual void GammaCorrection( const Types::DataItem gamma ) = 0;

  /// Function pointer type: double to double.
  typedef double (*FunctionTypeDouble)(const double);

  /// Function pointer type: float to float.
  typedef float (*FunctionTypeFloat)(const float);

  /** Apply real function to data. */ 	 
  virtual void ApplyFunction( Self::FunctionTypeFloat f ) 	 
  { 	 
    this->ApplyFunctionFloat( f ); 	 
  } 	 
  
  /** Apply real function to data. */ 	 
  virtual void ApplyFunction( Self::FunctionTypeDouble f ) 	 
  { 	 
    this->ApplyFunctionDouble( f ); 	 
  } 	 
	  	 
  /** Apply real function to data.
   */
  virtual void ApplyFunctionFloat( Self::FunctionTypeFloat f ) = 0;

  /** Apply real function to data.
   */
  virtual void ApplyFunctionDouble( Self::FunctionTypeDouble f ) = 0;

  /// Convert all values to absolute values.
  virtual void MakeAbsolute() = 0;

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
  virtual void Rescale( const Types::DataItem scale, const Types::DataItem offset, const Types::DataItem truncLo, const Types::DataItem truncHi = CMTK_ITEM_MAX ) = 0;

  /** Threshold data.
   * All values above upper threshold are set to upper thrershold. All values
   * below lower threshold are set to lower threshold.
   */
  virtual void Threshold( const Types::DataItemRange& range ) = 0;

  /** Threshold data.
   * All values outside the threshold range are set to the Padding (padding)
   * value.
   */
  virtual void ThresholdToPadding( const Types::DataItemRange& range ) = 0;

  /** Prune histogram to trim noise.
   * This function trims values on the upper and lower end of the value range by
   * thresholding in such a way that the trimmed number of samples on either side
   * amounts to the relative fraction of samples based on a target histogram bin
   * count.
   */
  virtual void PruneHistogram( const bool pruneHi, const bool pruneLo, const size_t numberOfBinsTarget, const size_t numberOfBinsInternal = 1024 );

  /** Binarize array values with given threshold.
   * All values greater than threshold (default: zero) are set to one, all
   * values smaller or equal are set to zero.
   */
  virtual void Binarize( const Types::DataItem threshold = 0 ) = 0;

  /** Clone an existing TypedArray.
   * Memory is allocated and the source array copied item by item. Data is
   * read by the source object's Get method and stored by the destination
   * object's Set method. 
   */
  Self::SmartPtr Clone() const 
  {
    return Self::SmartPtr( this->CloneVirtual() );
  }
  
  /// Default constructor.
  TypedArray () 
  {
    DataSize = 0;
    PaddingFlag = false;
    FreeArray = false;
    this->m_DataClass = DATACLASS_GREY;
  }
  
  /** Destructor.
   * Just has to be declared virtual.
   */
  virtual ~TypedArray () {}
  
  /** Free data pointer.
   * Derived classes may have to overload this method to take care of other
   * memory management libraries.
   */
  virtual void FreeData() = 0;

  /** Release pointer to data array.
   * The pointer must remain valid until object is destructed itself.
   */
  void ReleaseDataPointer () 
  {
    FreeArray = false;
  }
  
  /** Return the number of array elements.
   *@return The number of array elements
   */
  size_t GetDataSize () const { return DataSize; }

  /** Return the array size in bytes.
   */
  size_t GetDataSizeBytes () const
  { 
    return this->GetItemSize() * DataSize;
  }

  /** Return address of the data array.
   *@return A pointer to the real data array.
   */
  virtual void* GetDataPtr( const size_t offset = 0) = 0;

  /** Return address of the data array.
   *@return A pointer to the real data array.
   */
  virtual const void* GetDataPtr( const size_t offset = 0) const = 0;

  /** Get part of the stored array.
   *@return A pointer to the buffer given is returned.
   */
  virtual Types::DataItem* GetSubArray( Types::DataItem *const, const size_t, const size_t, const Types::DataItem = 0 ) const = 0;

  /** Allocate memory and get part of the stored array.
   *@return A pointer to the newly allocated buffer is returned.
   */
  virtual Types::DataItem* GetSubArray( const size_t, const size_t, const Types::DataItem = 0 ) const = 0;

  /** Return the flag for padding data.
   *@return If return value is zero, all data stored in the array is valid.
   * A non-zero value indicates, that there are locations with no data present.
   */
  bool GetPaddingFlag () const { return PaddingFlag; }

  /// Clear padding flag. This effectively turns padded pixels into actual data.
  void ClearPaddingFlag()
  {
    this->PaddingFlag = false;
  }

  /// Set the specified array element to no data present.
  virtual void SetPaddingAt ( const size_t ) = 0;

  /// Select the value to mark non-existent data.
  virtual void SetPaddingValue ( const Types::DataItem paddingData ) = 0;

  /// Select the value to mark non-existent data.
  virtual void SetPaddingPtr ( const void* paddingData ) = 0;

  /// Replace Padding data (padding) with given value.
  virtual void ReplacePaddingData ( const Types::DataItem value = 0 ) = 0;

  /// Check for padding data at given location.
  virtual bool IsPaddingAt( const size_t index ) const = 0;

  /// Check for padding data or zero at given location.
  virtual bool IsPaddingOrZeroAt( const size_t index ) const = 0;

  /// Calculate and return minimum and maximum data value.
  virtual const Types::DataItemRange GetRange() const = 0;

  /// Calculate entropy of distribution of values in this array.
  virtual double GetEntropy( const bool fractional = CMTK_HISTOGRAM_DISCRETE, const int numberOfBins = 128 ) const = 0;

  /** Calculate entropy of values in this array using existing histogram.
   * By using an already existing histogram, this function ensures that the
   * same numbers and sizes of bins are used for repeated calls of this 
   * function, even when the actual data in this array has changed. On the 
   * other hand, there is of course a risk that the data may have changed in
   * such a way that the old histogram arrangement does not fit it anymore.
   *
   * In the process of entropy computation, the histogram is also filled with
   * the data from this array. It is therefore available for subsequent
   * operations after returning from this function.
   */
  virtual double GetEntropy( Histogram<unsigned int>& histogram ) const = 0;

  /** Calculate entropy of values in this array using existing histogram.
   *\see TypedArray::GetEntropy
   */
  virtual double GetEntropy( Histogram<double>& histogram, const bool fractional = CMTK_HISTOGRAM_DISCRETE ) const = 0;

  /** Calculate entropy of values in this array using existing histogram and kernel for Parzen windowing.
   *\see TypedArray::GetEntropy
   */
  virtual double GetEntropy( Histogram<double>& histogram, const double* kernel, const size_t kernelRadius  ) const = 0;

  /** Calculate statistics.
   * Results will be both zero if there is not data in the array.
   *@return The number of valid (i.e., non-padding) values that constitute the
   * given results.
   */
  virtual size_t GetStatistics ( Types::DataItem& mean, Types::DataItem& variance ) const = 0;

  /** Compute approximate percentile value from histogram.
   */
  virtual Types::DataItem GetPercentile
  ( const Types::DataItem percentile, //!< The percentile to be computed. Value must be between 0 and 1.
    const size_t nBins = 256 //!< Number of histogram bins for percentile estimation.
    ) const;

  /** Compute list of approximate percentile values from histogram.
   * This function calls GetPercentile for each value in the given input vector and puts
   * all resulting values into the output vector in the same order. The main advantage of
   * using this function is that it is more efficient as a single histogram is created to
   * compute all percentiles.
   */
  virtual std::vector<Types::DataItem> GetPercentileList
  ( const std::vector<Types::DataItem>& percentileList, //!< The list of percentiles to be computed. Each value must be between 0 and 1.
    const size_t nBins = 256 //!< Number of histogram bins for percentile estimation.
    ) const;

  /** Get data histogram.
   *@return A histogram object filled with the relative frequencies of values 
   * in this array.
   */
  virtual Histogram<unsigned int>::SmartPtr GetHistogram( const unsigned int numberOfBins ) const = 0;

  /** Set data block to constant value.
   */
  virtual void BlockSet( const Types::DataItem value, const size_t fromOffset, const size_t toOffset ) = 0;

  /** Fill entire array with one value.
    */
  virtual void Fill( const Types::DataItem value )
  {
    this->BlockSet( value, 0, this->GetDataSize() );
  }

  /** Copy data block to other array.
   * This is really just a convenience wrapper for ConvertSubArray().
   */
  virtual void BlockCopy( Self& target, const size_t toOffset, const size_t fromOffset, const size_t blockLength ) const 
  {
    this->ConvertSubArray( target.GetDataPtr( toOffset ), target.GetType(), fromOffset, blockLength );
  }

  /** Exchange two data blocks.
   * Internally, the data is copied in chunks of up to 2KB (currently), using
   * an 'automatic' intermediate buffer of that size.
   */
  virtual void BlockSwap( const size_t fromOffset, const size_t toOffset, const size_t blockLength );

  /** Revert order of items in a given range.
   *\attention This implementation only works for item sizes up to 16 bytes
   * (128 bits); above that, the internal buffer size needs to be increased.
   */
  virtual void BlockReverse( const size_t fromOffset, const size_t blockLength );

  /** Apply function class to the values of this array.
   */
  virtual void ApplyFunctionObject( const TypedArrayFunction& f ) = 0;

protected:
  /** Scalar data type ID.
   * This field is not yet actually used anywhere. It merely serves a debugging
   * purpose since it allows easy identification of the type of this array from
   * inspecting its data.
   */
  ScalarDataType m_DataType;

  /// If not zero, the object should free the allocated array.
  bool FreeArray;

  /// The size of the data array, i.e. the number of items allocated.
  size_t DataSize;

  /** If true, PaddingFlag indicates there are items with no data present.
   */
  bool PaddingFlag;

  /// Allocate array for the given number of elements.
  virtual void Alloc ( const size_t datasize ) = 0;

  /// Virtual clone function: only to be calle by Clone().
  virtual Self* CloneVirtual() const = 0; 
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypedArray_h_included_
