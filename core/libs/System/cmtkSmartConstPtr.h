/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkSmartConstPtr_h_included_
#define __cmtkSmartConstPtr_h_included_

#include <cmtkconfig.h>

#ifndef NULL
#  define NULL 0
#endif

#include <algorithm>
#include <cassert>

#include <System/cmtkSafeCounter.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Smart pointer with reference counting.
 */
template<class T>
class SmartConstPointer
{
public:
  /// This class instance.
  typedef SmartConstPointer<T> Self;

  /// The underlying raw pointer type.
  typedef const T* PointerType;

  /// Null object.
  static Self Null;
  
  /// Get current reference counter value: use with caution!
  unsigned int GetReferenceCount() const 
  { 
    return this->m_ReferenceCount->Get(); 
  }

  /** Construct from dumb pointer.
   * Note that you MUST NEVER use this constructor more than once for each
   * dumb pointer, other than NULL!
   */
  explicit SmartConstPointer( T *const object = NULL ) 
    : m_ReferenceCount( new SafeCounter( 1 ) )
  { 
    this->m_Object.ptrConst = object;
  }

  /** Copy constructor template.
   * Increment reference counter in the process.
   */
  template<class T2>
  SmartConstPointer( const SmartConstPointer<T2>& ptr ) 
    : m_ReferenceCount( ptr.m_ReferenceCount )
  {
    this->m_Object.ptrConst = ptr.m_Object.ptrConst;
    this->m_ReferenceCount->Increment();
  }
  
  /** Copy constructor to prevent compiler-generated copy constructor.
   * Increment reference counter in the process.
   */
  SmartConstPointer( const Self& ptr ) 
    : m_ReferenceCount( ptr.m_ReferenceCount )
  {
    this->m_Object.ptrConst = ptr.m_Object.ptrConst;
    this->m_ReferenceCount->Increment();
  }
  
  /// Destruct: decrease reference pointer and free dumb pointer if necessary.
  ~SmartConstPointer() 
  {
    assert( this->m_ReferenceCount != NULL ); // we may have m_Object=NULL, but m_ReferenceCount should never be!
    if ( ! this->m_ReferenceCount->Decrement() ) 
      {
      delete this->m_ReferenceCount;
#ifdef DEBUG
      this->m_ReferenceCount = NULL;
#endif
      if ( this->m_Object.ptrConst ) 
	{
	delete this->m_Object.ptrConst;
#ifdef DEBUG
	this->m_Object.ptrConst = NULL;
#endif
	}
      }
  }

  /// De-referencing operator (returns constant object).
  const T& operator*() const { return *this->m_Object.ptrConst; }

  /// De-referencing operator (returns constant object pointer).
  const T* operator->() const { return this->m_Object.ptrConst; }

  /// Implicit conversion to constant pointer.
  operator const T*() const { return m_Object.ptrConst; }
  
  /// Explicit conversion to constant pointer.
  const T* GetConstPtr() const { return this->m_Object.ptrConst; }

  /** Release control of this pointer.
   *\note This is a dangerous function. Be sure you know what you are doing!
   */
  const T* ReleasePtr() 
  { 
    const T* object = this->m_Object.ptrConst; 
    this->m_Object.ptrConst = NULL; 
    return object; 
  }

  /** Assignment operator.
   * This is implemented using the std::swap function.
   *\warning The "other" parameter HAS TO USE CALL BY VALUE for this function to work,
   *  because we are not creating an explicit copy of the original object before 
   *  calling Swap() (see Effective C++, 3rd, Item 11, p.56).
   *\warning Open question: given that we pass the parameter by value, not reference, does
   *  this really prevent the definition of a compiler-generated assignment (which passes
   *  parameter by reference)?
   */
  const Self& operator= ( const Self other ) const
  {
    using std::swap;
    swap( this->m_ReferenceCount, other.m_ReferenceCount );
    swap( this->m_Object.ptrConst, other.m_Object.ptrConst );

    return *this;
  }

  /** Equality operator using pointer.
   */
  bool operator== ( const Self& other ) const 
  {
    return (this->m_Object.ptrConst == other.m_Object.ptrConst);
  }
  
  /** Inequality operator.
   */
  bool operator!= ( const Self& other ) const 
  {
    return (this->m_Object.ptrConst != other.m_Object.ptrConst);
  }
  
  /** "Smaller than" operator (necessay for storing objects in std::map).
   */
  bool operator< ( const Self& other ) const 
  {
    return ( this->m_Object.ptrConst < other.m_Object.ptrConst );
  }
  
  ///Dynamic cast between smart pointer types.
  template<class T2> 
  static Self DynamicCastFrom( const T2& from_P )
  {
    return Self( dynamic_cast<typename Self::PointerType>( from_P.GetConstPtr() ), from_P.m_ReferenceCount );
  }

protected:
  /// Pointer to detached reference counter for this object.
  mutable SafeCounter* m_ReferenceCount;
  
  /// Union for const and non-const pointer for reference-counted object.
  mutable union
  {
    /// Pointer to const.
    const T* ptrConst;
    /// Pointer to non-const.
    T* ptr;
  } m_Object;
  
  /** Construct from dumb pointer and existing reference counter.
   * The reference counter is increased in the process.
   */
  SmartConstPointer( const T *const object, SafeCounter *const counter ) 
  {
    this->m_Object.ptrConst = object;
    this->m_ReferenceCount = counter;
    this->m_ReferenceCount->Increment();
  }
  
  /// Make all template instances friends for easy type casting.
  template<class T2> friend class SmartConstPointer;
};

template<typename T> SmartConstPointer<T> SmartConstPointer<T>::Null;

//@}

} // namespace cmtk

#endif // define __cmtkSmartConstPtr_h_included_
