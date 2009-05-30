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

#ifndef __cmtkVector_h_included_
#define __cmtkVector_h_included_

#include <cmtkconfig.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>

#include <algorithm>

#include <cmtkMathUtil.h>
#include <cmtkException.h>
#include <cmtkTypes.h>
#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Error message used when dimensions of two vectors do not match.
extern const char* DimensionMismatchError;

/** Numerical vector class.
 *@author Torsten Rohlfing
 */
template<class T>
class Vector
{
public:
  /// Vector dimension.
  size_t Dim;

  /// Vector elements.
  T *Elements;

private:
  /// Flag for memory deallocation of value array.
  bool FreeElements;

public:
  /// Smart pointer to igsFloatVector.
  typedef SmartPointer< Vector<T> > SmartPtr;

  /**@name Constructors */
  //@{
  /// Create constant (zero-)vector.
  Vector ( const size_t dim = 0, const T value = 0 ) 
  {
    Dim = dim;
    if ( Dim ) 
      {
      Elements = Memory::AllocateArray<T>( Dim );
      FreeElements = true;
      if ( value==0 )
	memset( Elements, 0, Dim * sizeof(T) );
      else
	for ( size_t i=0; i<Dim; ++i )
	  Elements[i]=value;
      } 
    else
      {
      Elements = NULL;
      FreeElements = false;
      }
  }

  /** Create vector from existing array.
   */
  Vector ( const size_t dim, T *const elements, const bool freeElements = true ) 
  {
    Dim = dim;
    Elements = elements;
    FreeElements = freeElements;
  }
  
  /// Create vector from other vector (also subvector).
  Vector ( const Vector& other, const size_t len = 0, const size_t from = 0 ) 
  {
    if ( len )
      Dim = std::min( len, other.Dim - from );
    else
      Dim = other.Dim - from;
    
    Elements = Memory::AllocateArray<T>( Dim );
    FreeElements = true;
    memcpy( Elements, other.Elements + from, Dim * sizeof(T) );
  }
  //@}

  /// Clone (sub)vector.
  Vector* Clone( const size_t len = 0, const size_t from = 0 ) const
  { 
    return new Vector( *this, len, from ); 
  }
  
  /// Destructor.
  ~Vector () 
  {
    if ( Elements && FreeElements ) 
      {
      delete[] Elements;
      }
  }
  
  /** Check vector dimension.
   * This object's dimension as a vector is compared to a given value. An
   * exception is thrown if both values are different.
   */
  void CheckDim ( const size_t dim ) const 
  {
    if ( dim != Dim )
      throw Exception( DimensionMismatchError, this );
  }

  /** Test vector dimension.
   * This object's dimension as a vector is compared to a given value. Other
   * than CheckDim, this method does not throw an excpetion on a mismatch.
   */
  bool CompareDim ( const size_t dim ) const 
  {
    return ( dim == Dim );
  }

  /** Set vector dimension.
   * If the current vector dimension is not equal to the requested dimension,
   * the elements array is deleted and a new one is allocated. In any case,
   * there is no guarantee that the data stored in the vector before this call
   * remains unchanged. This is even true for initial elements.
   *@param dim The number of elements to be stored in this vector after
   * returning from this function.
   *@param zero If this parameter is true, all vector elements are set to
   * the zero value in their respective data type.
   *@return A reference to this object after changing the dimension.
   */
  Vector& SetDim ( const size_t dim, const bool zero = true ) 
  {
    if ( Dim != dim ) 
      {
      if ( Elements ) 
	{
	delete[] Elements;
	}

      Dim = dim;
      
      if ( Dim ) 
	{
	Elements = Memory::AllocateArray<T>( Dim );
	} 
      else
	Elements = NULL;
      }
    
    if ( zero && Dim ) 
      {
      memset( Elements, 0, Dim * sizeof(T) );
      }
    
    return *this;
  }
  
  /** Adjust vector dimension.
   * Unlike SetDim(), this function preserves the values of elements in the
   * vector if they are still in the valid index range after size adjustment.
   *@param dim The number of elements to be stored in this vector after
   * returning from this function.
   *@param zero If this parameter is true, all new vector elements are set to
   * the zero value in their respective data type.
   *@return A reference to this object after changing the dimension.
   */
  Vector& AdjustDimension( const size_t dim, const bool zero = true ) 
  {
    // If old and new size are the same, there is nothing to do.
    if ( Dim != dim ) 
      {
      T* newElements = Memory::AllocateArray<T>( dim );
      // copy common elements
      memcpy( newElements, this->Elements, sizeof(T) * std::min( dim, Dim ) );

      // reset new elements if so desired
      if ( zero && (dim > Dim) )
	{
	memset( newElements + Dim, 0, sizeof(T) * (dim-Dim) );
	}

      // new set new array.
      this->Dim = dim;
      this->Elements = newElements;
      } 
    
    return *this;
  }
  
  /// Vector assignment.
  Vector& operator = ( const Vector& other ) 
  {
    if ( Dim != other.Dim ) {
    if (Elements) 
      {
      delete[] Elements;
      Elements = NULL;
      }
    
    Dim = other.Dim;
    }
    
    if ( Elements == NULL ) 
      {
      Elements = Memory::AllocateArray<T>( Dim );
      }
    
    memcpy( Elements, other.Elements, Dim * sizeof(T) );
    return *this; 
  }

  /** Copy another vector to given offset.
   *@param other Vector from which the specified elements are copied.
   *@param offs Destination offset. Copying starts at this position in this
   * instance.
   *@param len Number of elements to be copied. If zero, all elements are 
   * copied until the end of one of the vectors is reached.
   */
  void CopyToOffset( const Vector& other, const size_t offs = 0, size_t len = 0 )
  {
    if ( ! len ) 
      len = std::min( this->Dim - offs, other.Dim );
    for ( size_t idx=0; idx<len; ++idx )
      Elements[offs+idx] = other.Elements[idx];
  }

  /// Test for vector equality.
  int operator== ( const Vector& other ) const 
  {
    if ( Dim != other.Dim )
      return 0;
    
    for ( size_t i=0; i<Dim; ++i )
      if ( Elements[i] != other.Elements[i] )
	return 0;
    
    return 1;
  }
  
  /// Test for vector inequality.
  int operator!= ( const Vector& other ) const 
  {
    return !(*this == other );
  }
  
  /// Calculate Euclid's vector norm.
  T EuclidNorm () const 
  { 
    T Result = 0;
    
#pragma omp parallel for if (Dim>1e4) reduction(+:Result)
    for ( size_t i=0; i<this->Dim; ++i ) 
      {
      const T e = Elements[i];
      Result+=e*e;
      }
    
    return sqrt(Result);
  }

  /// Calculate maximum vector norm.
  T MaxNorm () const 
  { 
    T Result = 0;
    
    for ( size_t i=0; i<Dim; ++i ) 
      {
      Result = std::max<T>( Result, fabs( Elements[i] ) );
      }
    
    return Result;
  }

  /// Set all vector elements to zero.
  void Clear() 
  { 
    memset( Elements, 0, Dim * sizeof( *Elements ) ); 
  }

  /// Set all vector elements to constant value.
  void SetAll( const T value )
  {
#pragma omp parallel for if (Dim>1e5)
    for ( size_t i=0; i < this->Dim; ++i ) 
      this->Elements[i] = value;
  }

  /// Test if vector is all zero.
  int IsNull () const 
  { 
    for ( size_t i=0; i<Dim; ++i )
      if ( Elements[i] != 0 )
	return 0;
    
    return 1;
  }  
    
  /// Get vector element by coordinate index.
  T& operator [] ( const size_t index ) 
  {
    return Elements[index];
  }

  /// Get constant vector element by coordinate index.
  const T& operator [] ( const size_t index ) const 
  {
    return Elements[index];
  }

  /// Increment vector by another.
  Vector<T>& operator+= ( const Vector<T>& delta ) 
  {
    assert( Dim == delta.Dim );

#pragma omp parallel for if (Dim>1e4)
    for ( size_t i=0; i<this->Dim; ++i )
      Elements[i] += delta.Elements[i];
    
    return *this;
  }

  /// Increment vector by another, posibly with a scalar weight.
  Vector<T>& Add( const Vector<T>& delta, const T weight = 1 ) 
  {
    assert( Dim == delta.Dim );

#pragma omp parallel for if (Dim>1e4)
    for ( size_t i=0; i< this->Dim; ++i )
      Elements[i] += weight * delta.Elements[i];
    
    return *this;
  }

  /// Set partial vector by copying from another.
  Vector<T>& SetPartial 
  ( const Vector<T>& other, const size_t offs = 0 ) 
  {
    for ( size_t i=0; (i<other.Dim) && (i+offs<Dim); ++i )
      Elements[offs+i] = other.Elements[i];
    
    return *this;
  }

  /// Increment partial vector by another.
  Vector<T>& AddPartial 
  ( const Vector<T>& delta, const size_t offs = 0 ) 
  {
    for ( size_t i=0; (i<delta.Dim) && (i+offs<Dim); ++i )
      Elements[offs+i]+=delta.Elements[i];
    
    return *this;
  }

  /// Decrement vector by another.
  Vector<T>& operator-= ( const Vector<T>& delta ) 
  {
    assert( Dim == delta.Dim );
    
#pragma omp parallel for if (Dim>1e4)
    for ( size_t i=0; i < this->Dim; ++i )
      Elements[i] -= delta.Elements[i];
    
    return *this;
  }

  /// Subtract partial vector from another.
  Vector<T>& SubtractPartial 
  ( const Vector<T>& delta, const size_t offs = 0 ) 
  {
    for ( size_t i=0; (i<delta.Dim) && (i+offs<Dim); ++i )
      this->Elements[i+offs]-=delta.Elements[i];
    
    return *this;
  }

  /// Multiply by a scalar.
  Vector<T>& operator*= ( const T a ) 
  {
#pragma omp parallel for if (Dim>1e4)
    for ( size_t i=0; i<this->Dim; ++i )
      this->Elements[i] *= a;
    
    return *this;
  }

  void Print ( FILE *const fp = NULL, const char* format = " %f" ) const 
  {
    if ( fp ) 
      {
      for ( size_t idx=0; idx < Dim; ++idx )
	fprintf( fp, format, (float) Elements[idx] );
      fputs( "\n", fp );
      } 
    else
      {
      for ( size_t idx=0; idx < Dim; ++idx )
	printf( format, (float) Elements[idx] );
      puts( "" );
      }
  }
  
  /** Sort values in the vector.
   * Using the two parameters, from and len, this function can be used to sort
   * only a subrange of values in this vector. In particular, it can be used to 
   * sort the first len elements if from == 0.
   *\param from Index of first element in the range to sort.
   *\param len Number of elements to be sorted.
   */
  void Sort( const size_t from = 0, const size_t len = 0 ) 
  {
    T *ptr = Elements+from;
    if ( len )
      qsort( ptr, len, sizeof( T ), Vector<T>::Compare );
    else
      qsort( ptr, Dim-from, sizeof( T ), Vector<T>::Compare );
  }
  
private:
  /// Compare two vector elements; this is needed for sorting.
  static int Compare( const void* a, const void* b ) 
  {
    const T *Ta = (const T *) a;
    const T *Tb = (const T *) b;
    return (*Ta > *Tb) - (*Ta < *Tb);
  }
};

/** Shortcut definition.
 * This typedef defines a name for the frequently used vectors over the 
 * Types::Coordinate type. This is used for all kinds of parameters vectors etc.
 */
typedef Vector<Types::Coordinate> CoordinateVector;

/** Shortcut definition.
 * This typedef defines a name for the frequently used vectors over the 
 * float type.
 */
typedef Vector<float> FloatVector;

//@}

} // namespace cmtk

#include <cmtkVector.txx>

#endif // #ifndef __cmtkVector_h_included_
