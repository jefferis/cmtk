/*
//
//  Copyright 2010-2012 SRI International
//
//  Copyright 2010 Torsten Rohlfing
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

#ifndef __cmtkFixedArray_h_included_
#define __cmtkFixedArray_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkDataTypeTraits.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

#include <algorithm>
#include <iostream>

namespace
cmtk
{
/// Class for fixed-size n-dimensional array.
template<size_t NDIM,typename T=int>
class FixedArray
{
public:
  /// This class.
  typedef FixedArray<NDIM,T> Self;

  /// Type of the stored values.
  typedef T ValueType;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer-to-const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Return vector size.
  size_t Size() const { return NDIM; }

  /// Zero operator.
  static const Self Zero()
  {
    Self v;
    std::fill( v.begin(), v.end(), DataTypeTraits<T>::Zero() );
    return v;
  }

  /** Default constructor. 
    * note for improved efficiency, this does NOT initialize the data array.
    */
  FixedArray() {}

  /// Initialization constructor.
  explicit FixedArray( const T& initValue )
  {
    std::fill( this->begin(), this->end(), initValue );
  }

  /// Type conversion constructor template.
  template<class T2>
  FixedArray( const FixedArray<NDIM,T2>& rhs )
  {
    for ( size_t i = 0; i < NDIM; ++i )
      this->m_Data[i] = static_cast<T>( rhs[i] );
  }

  /// Make array from const pointer.
  template<class T2> static Self FromPointer( const T2 *const ptr ) 
  { 
    Self v;
    for ( size_t i = 0; i < NDIM; ++i )
      v[i] = ptr[i];

    return v;
  }

  /// Get element reference.
  T& operator[]( const size_t i )
  {
    return this->m_Data[i];
  }

  /// Get const element reference.
  const T& operator[]( const size_t i ) const
  {
    return this->m_Data[i];
  }

  /// Equality operator.
  bool operator==( const Self& rhs ) const
  {
    for ( size_t i = 0; i<NDIM; ++i )
      if ( this->m_Data[i] != rhs.m_Data[i] )
	return false;
    return true;
  }

  /// Inequality operator.
  bool operator!=( const Self& rhs ) const
  {
    return !this->operator==( rhs );
  }

  /// Pointer to first array element.
  T* begin()
  {
    return this->m_Data;
  }

  /// Pointer behind last array element.
  T* end()
  {
    return this->m_Data+NDIM;
  }

  /// Pointer to first array element.
  const T* begin() const
  {
    return this->m_Data;
  }

  /// Pointer behind last array element.
  const T* end() const
  {
    return this->m_Data+NDIM;
  }

protected:
  /// The actual index array.
  T m_Data[NDIM];
};

/// Stream output operator.
template<size_t NDIM,typename T>
std::ostream& operator<<( std::ostream& stream, const FixedArray<NDIM,T>& index )
{
  for ( size_t i = 0; i < NDIM; ++i )
    stream << index[i] << " ";

  return stream;
}

/// Stream input operator.
template<size_t NDIM,typename T>
std::istream& operator>>( std::istream& stream, FixedArray<NDIM,T>& index )
{
  for ( size_t i = 0; i < NDIM; ++i )
    stream >> index[i];
  return stream;
}

} // namespace cmtk

#endif // #ifndef __cmtkFixedArray_h_included_
