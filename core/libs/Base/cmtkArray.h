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

#ifndef __cmtkArray_h_included_
#define __cmtkArray_h_included_

#include <cmtkconfig.h>

#include <cmtkMemory.h>

#include <stddef.h>
#include <stdlib.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/// One-dimensional array template.
template<class T>
class Array
{
public:
  /// Constructor: allocate array.
  Array( const size_t dim ) 
  {
    this->m_Size = dim;
    this->m_Array = Memory::AllocateArray<T>( dim );
  }

  /// Destructor: free allocated array.
  ~Array() { delete[] this->m_Array; }

  /// Get dimension of this array.
  size_t GetDim() const { return this->m_Size; }

  /// Return pointer to constant array.
  operator const T*() const { return this->m_Array; }

  /// Return pointer to array.
  operator T*() { return this->m_Array; }

  /// Reset all values to zero.
  void SetAllToZero() 
  {
    memset( this->m_Array, 0, this->m_Size * sizeof( T ) );
  }

  /// Set all values.
  void SetAll( const T value ) 
  {
    for ( size_t i = 0; i < this->m_Size; ++i ) 
      {
      this->m_Array[i] = value;
      }
  }

  /// Copy values from another array.
  Array<T>& operator= ( const Array<T>& other ) 
  {
    memcpy( this->m_Array, other.m_Array, this->m_Size * sizeof( T ) );
    return *this;
  }

protected:
  /// The actual array.
  T* m_Array;

  /// The number of elements allocated.
  size_t m_Size;
};

} // namespace

#endif // #ifndef __cmtkArray_h_included_
