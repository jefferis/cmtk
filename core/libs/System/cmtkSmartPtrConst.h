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

#ifndef __cmtkSmartPtrConst_h_included_
#define __cmtkSmartPtrConst_h_included_

#include <cmtkconfig.h>

#ifndef NULL
#  define NULL 0
#endif

#include <assert.h>

#include <cmtkSafeCounter.h>
#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Smart pointer with reference counting.
 */
template<class T>
class SmartPointerConst
{
private:
  /// Pointer to detached reference counter for this object.
  mutable SafeCounter* ReferenceCount;
  
  /// Pointer for reference-counted object.
  mutable T* Object;
  
  /** Construct from dumb pointer and existing reference counter.
   * The reference counter is increased in the process.
   */
  SmartPointerConst( T *const object, SafeCounter *const counter ) 
  {
    Object = object;
    ReferenceCount = counter;
    ReferenceCount->Increment();
  }
  
public:
  /// This class instance.
  typedef SmartPointerConst<T> Self;

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
  explicit SmartPointerConst( T *const object = NULL ) 
  { 
    ReferenceCount = new SafeCounter( 1 );
    Object = object;
  }

  /** Construct from const smart pointer reference.
   * Increment reference counter in the process.
   */
  SmartPointerConst( const Self& ptr ) 
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
  
  /** Construct from non-const smart pointer reference.
   * Increment reference counter in the process.
   */
  SmartPointerConst( const SmartPointer<T>& ptr ) 
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
  ~SmartPointerConst() 
  {
    assert( ReferenceCount != NULL ); // we may have Object=NULL, but ReferenceCount should never be!
    if ( ! ReferenceCount->Decrement() ) 
      {
      delete ReferenceCount;
      if ( Object ) 
	{
	delete Object;
	}
      }
  }

  /// De-referencing operator (returns constant object).
  const T& operator*() const { return *Object; }

  /// De-referencing operator (returns volatile object).
  const volatile T& operator*() const volatile { return *Object; }

  /// De-referencing operator (returns constant object pointer).
  const T* operator->() const { return Object; }

  /// De-referencing operator (returns constant object pointer).
  const volatile T* operator->() const volatile { return Object; }

  /// Implicit conversion to constant pointer.
  operator const T*() const { return Object; }
  
  /// Explicit conversion to bool (validity).
  bool IsNull() const { return Object == NULL; }
  
  /// Explicit conversion to constant pointer.
  const T* GetPtr() const { return Object; }

  /// Explicit conversion to constant and volatile pointer.
  const volatile T* GetPtr() const volatile { return Object; }

  /** Release control of this pointer.
   *\note This is a dangerous function. Be sure you know what you are doing!
   */
  const T* ReleasePtr() 
  { 
    const T* object = Object; 
    Object = NULL; 
    return object; 
  }

  /** Assignment operator.
   * Implemented using Swap() operator.
   */
  Self& operator= ( const Self& other ) 
  {
    Self temp( other );
    this->Swap( temp );
    return *this;
  }

  /** Assignment operator from non-const smart pointer.
   * Implemented using Swap() operator.
   */
  Self& operator= ( const SmartPointer<T>& other ) 
  {
    Self temp( other );
    this->Swap( temp );
    return *this;
  }

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
  template<class T2> operator SmartPointerConst<T2>() 
  { 
    return SmartPointerConst<T2>( Object, ReferenceCount );
  }
  
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
  void Swap( Self& other ) 
  {
    this->SwapPrimitive( ReferenceCount, other.ReferenceCount );
    this->SwapPrimitive( Object, other.Object );
  }
  
  /// Helper function that swaps two primitive objects (or pointers).
  template<class TT> void SwapPrimitive( TT& a, TT& b ) 
  {
    TT temp = a;
    a = b;
    b = temp;
  }

  /// Make all template instances friends for easy type casting.
  template<class T2> friend class SmartPointerConst;
};

template<typename T> SmartPointerConst<T> SmartPointerConst<T>::Null;

//@}

} // namespace cmtk

#endif // define __cmtkSmartPtrConst_h_included_
