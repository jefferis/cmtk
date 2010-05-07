/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkFixedVector_h_included_
#define __cmtkFixedVector_h_included_

#include <cmtkconfig.h>

#include <fstream>

namespace
cmtk
{
/// Class for n-dimensional image index.
template<size_t NDIM,typename T=int>
class FixedVector
{
public:
  /// This class.
  typedef FixedVector<NDIM,T> Self;

  /// Type of the stored values.
  typedef T ValueType;

  /// Default constructor.
  FixedVector() {}

  /// Constructor from const int array.
  FixedVector( const T (&indexArray)[NDIM] )
  {
    for ( size_t i = 0; i < NDIM; ++i )
      this->m_FixedVector[i] = indexArray[i];
  }

  /// Get element reference.
  T& operator[]( const size_t i )
  {
    return this->m_FixedVector[i];
  }

  /// Get const element reference.
  const T& operator[]( const size_t i ) const
  {
    return this->m_FixedVector[i];
  }

  /// Addition operator.
  const Self operator+( const Self& rhs ) const
  {
    return Self( *this ) += rhs;
  }
  
  /// In-place addition operator.
  Self& operator+=( const Self& rhs )
  {
    for ( size_t i = 0; i<NDIM; ++i )
      this->m_FixedVector[i] += rhs.m_FixedVector[i];
    return *this;
  }
  
  /// Subtraction operator.
  const Self operator-( const Self& rhs ) const
  {
    return Self( *this ) -= rhs;
  }
  
  /// In-place subtraction operator.
  Self& operator-=( const Self& rhs )
  {
    for ( size_t i = 0; i<NDIM; ++i )
      this->m_FixedVector[i] -= rhs.m_FixedVector[i];
    return *this;
  }

  /// Equality operator.
  bool operator==( const Self& rhs ) const
  {
    for ( size_t i = 0; i<NDIM; ++i )
      if ( this->m_FixedVector[i] != rhs.m_FixedVector[i] )
	return false;
    return true;
  }

  /// Pointer to first array element.
  T* begin()
  {
    return m_FixedVector;
  }

  /// Pointer behind last array element.
  T* end()
  {
    return m_FixedVector+NDIM;
  }

  /// Pointer to first array element.
  const T* begin() const
  {
    return m_FixedVector;
  }

  /// Pointer behind last array element.
  const T* end() const
  {
    return m_FixedVector+NDIM;
  }

private:
  /// The actual index array.
  T m_FixedVector[NDIM];
};

/// Stream input operator.
template<size_t NDIM,typename T>
std::ofstream& operator<<( std::ofstream& stream, const FixedVector<NDIM,T>& index )
{
  for ( size_t i = 0; i < NDIM; ++i )
    stream << index[i];
  return stream;
}

/// Stream output operator.
template<size_t NDIM,typename T>
std::ifstream& operator>>( std::ifstream& stream, FixedVector<NDIM,T>& index )
{
  for ( size_t i = 0; i < NDIM; ++i )
    stream >> index[i];
  return stream;
}

} // namespace cmtk

#endif // #ifndef __cmtkFixedVector_h_included_
