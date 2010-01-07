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

#ifndef __cmtkSmartPtr_h_included_
#define __cmtkSmartPtr_h_included_

#include <cmtkconfig.h>

#ifndef NULL
#  define NULL 0
#endif

#ifdef DEBUG
#  include <stdio.h>
#endif
#include <assert.h>

#include <cmtkSafeCounter.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Smart pointer with reference counting.
 */
template<class T>
class SmartPointer
{
private:
  /// Pointer to detached reference counter for this object.
  mutable SafeCounter* ReferenceCount;
  
  /// Pointer for reference-counted object.
  mutable T* Object;
  
  /** Construct from dumb pointer and existing reference counter.
   * The reference counter is increased in the process.
   */
  SmartPointer( T *const object, SafeCounter *const counter ) 
  {
    Object = object;
    ReferenceCount = counter;
    ReferenceCount->Increment();
  }
  
public:
  /// This class instance.
  typedef SmartPointer<T> Self;

  /// The underlying raw pointer type.
  typedef T* PointerType;

  /// Null object.
  static Self Null;
  
  /// Get current reference counter value: use with caution!
  unsigned int GetReferenceCount() const { return ReferenceCount->Get(); }

  /** Construct from dumb pointer.
   * Note that you MUST NEVER use this constructor more than once for each
   * dumb pointer, other than NULL!
   */
  explicit SmartPointer( T *const object = NULL ) 
  { 
    ReferenceCount = new SafeCounter( 1 );
    Object = object;
  }

  /** Construct from other smart pointer reference.
   * Increment reference counter in the process.
   */
  SmartPointer( const Self& ptr ) 
  {
    if ( &ptr ) 
      {
      ReferenceCount = ptr.ReferenceCount;
      ReferenceCount->Increment();
      Object = ptr.Object;
      } 
    else 
      {
      ReferenceCount = new SafeCounter( 1 );
      Object = NULL;
      }
  }
  
  /// Destruct: decrease reference pointer and free dumb pointer if necessary.
  ~SmartPointer() 
  {
    assert( ReferenceCount != NULL ); // we may have Object=NULL, but ReferenceCount should never be!
    if ( ! ReferenceCount->Decrement() ) 
      {
      delete ReferenceCount;
#ifdef DEBUG
      ReferenceCount = NULL;
#endif
      if ( Object ) 
	{
	delete Object;
#ifdef DEBUG
	Object = NULL;
#endif
	}
      }
  }

  /// De-referencing operator (returns non-constant object).
  T& operator*() { return *Object; }

  /// De-referencing operator (returns constant object).
  const T& operator*() const { return *Object; }

  /// De-referencing operator (returns non-constant object pointer).
  T* operator->() { return Object; }

  /// De-referencing operator (returns constant object pointer).
  const T* operator->() const { return Object; }

  /// Implicit conversion to constant pointer.
  operator const T*() const { return Object; }
  
  /// Explicit conversion to bool (validity).
  bool IsNull() const { return Object == NULL; }
  
  /// Explicit conversion to non-constant pointer.
  T* GetPtr() { return Object; }

  /// Explicit conversion to constant pointer.
  const T* GetPtr() const { return Object; }

  /** Release control of this pointer.
   *\note This is a dangerous function. Be sure you know what you are doing!
   */
  T* ReleasePtr() 
  { 
    T* object = Object; 
    Object = NULL; 
    return object; 
  }

  /** Assignment operator.
   * Implemented using Swap() operator.
   */
  const Self& operator= ( const Self& other ) const
  {
    const Self temp( other );
    this->Swap( temp );
    return *this;
  }

#ifndef __SUNPRO_CC
#ifndef _MSC_VER
  /** Assignment operator.
   * Implemented using Swap() operator.
   */
  Self& operator= ( Self& other )
  {
    Self temp( other );
    this->Swap( temp );
    return *this;
  }
#endif
#endif

  /** Equality operator using pointer.
   */
  bool operator== ( const Self& other ) const 
  {
    return (this->Object == other.Object);
  }
  
  /** Inequality operator.
   */
  bool operator!= ( const Self& other ) const 
  {
    return (this->Object != other.Object);
  }
  
  /** "Smaller than" operator (necessay for storing objects in std::map).
   */
  bool operator< ( const Self& other ) const 
  {
    return ( this->Object < other.Object );
  }
  
  /// Implicit cast conversion operator.
  template<class T2> operator SmartPointer<T2>() 
  { 
    return SmartPointer<T2>( Object, ReferenceCount );
  }
  
#ifndef _MSC_VER
  /// Implicit cast conversion operator.
  template<class T2> operator const SmartPointer<T2>() const
  { 
    return SmartPointer<T2>( Object, ReferenceCount );
  }
#endif
  
  ///Dynamic cast between smart pointer types.
  template<class T2> 
  static Self DynamicCastFrom( const T2& from_P )
  {
    Self result( dynamic_cast<typename Self::PointerType>( from_P.Object ), from_P.ReferenceCount );
    return result;
  }

private:
  /** Swap two smart pointers.
   * This function is used for reference-safe assignment (i.e., replacing) of
   * smart pointer objects.
   */
  void Swap( const Self& other ) const
  {
    Self::SwapPrimitive( ReferenceCount, other.ReferenceCount );
    Self::SwapPrimitive( Object, other.Object );
  }
  
  /** Swap two smart pointers.
   * This function is used for reference-safe assignment (i.e., replacing) of
   * smart pointer objects.
   */
  void Swap( Self& other )
  {
    Self::SwapPrimitive( ReferenceCount, other.ReferenceCount );
    Self::SwapPrimitive( Object, other.Object );
  }
  
  /// Helper function that swaps two primitive objects (or pointers).
  template<class TT> static void SwapPrimitive( TT& a, TT& b )
  {
    TT temp = a;
    a = b;
    b = temp;
  }

  /// Make all template instances friends for easy type casting.
  template<class T2> friend class SmartPointer;
};

template<typename T> SmartPointer<T> SmartPointer<T>::Null;

//@}

} // namespace cmtk

#endif // define __cmtkSmartPtr_h_included_
