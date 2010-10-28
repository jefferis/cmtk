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

#ifndef __cmtkSmartPtr_h_included_
#define __cmtkSmartPtr_h_included_

#include <cmtkconfig.h>

#ifndef NULL
#  define NULL 0
#endif

#include <algorithm>
#include <cassert>

#include <System/cmtkSmartConstPtr.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Smart pointer with reference counting.
 */
template<class T>
class SmartPointer :
  /// Inherit from smart pointer-to-const.
  public SmartConstPointer<T>
{
public:
  /// This class instance.
  typedef SmartPointer<T> Self;

  /// This class instance.
  typedef SmartConstPointer<T> Superclass;

  /// Use "const" inherited member functions.
  using Superclass::operator*;
  using Superclass::operator->;
  using Superclass::ReleasePtr;

  /// The underlying raw pointer type.
  typedef T* PointerType;

  /// Null object.
  static Self Null;
  
  /** Default constructor.
   */
  SmartPointer() {}
  
  /** Construct from dumb pointer.
   * Note that you MUST NEVER use this constructor more than once for each
   * dumb pointer, other than NULL!
   */
  explicit SmartPointer( T *const object ) : SmartConstPointer<T>( object ) {}

  /** Copy constructor template.
   * Increment reference counter in the process.
   */
  template<class T2>
  SmartPointer( const SmartPointer<T2>& ptr ) : SmartConstPointer<T>( ptr ) {}
  
  /** Copy constructor to prevent compiler-generated copy constructor.
   * Increment reference counter in the process.
   */
  SmartPointer( const Self& ptr ) : SmartConstPointer<T>( ptr ) {}
  
  /// De-referencing operator (returns non-constant object).
  T& operator*() { return *this->m_Object.ptr; }

  /// De-referencing operator (returns non-constant object pointer).
  T* operator->() { return this->m_Object.ptr; }

  /// Explicit conversion to non-constant pointer.
  T* GetPtr() const { return this->m_Object.ptr; }

  /** Release control of this pointer.
   *\note This is a dangerous function. Be sure you know what you are doing!
   */
  T* ReleasePtr() 
  { 
    T* object = this->m_Object.ptr; 
    this->m_Object.ptr = NULL; 
    return object; 
  }

  ///Dynamic cast between smart pointer types.
  template<class T2> 
  static Self DynamicCastFrom( const T2& from_P )
  {
    return Self( dynamic_cast<typename Self::PointerType>( from_P.GetPtr() ), from_P.m_ReferenceCount );
  }

private:
  /** Construct from dumb pointer and existing reference counter.
   * The reference counter is increased in the process.
   */
  SmartPointer( T *const object, SafeCounter *const counter ) : Superclass( object, counter ) {}
  
  /// Make all template instances friends for easy type casting.
  template<class T2> friend class SmartPointer;
};

template<typename T> SmartPointer<T> SmartPointer<T>::Null;

//@}

} // namespace cmtk

#endif // define __cmtkSmartPtr_h_included_
