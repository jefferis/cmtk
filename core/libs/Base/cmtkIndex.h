/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkIndex_h_included_
#define __cmtkIndex_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{
/// Class for n-dimensional image index.
template<size_t NDIM,typename T=int>
class Index
{
public:
  /// This class.
  typedef Index<NDIM,T> Self;

  /// Type of the stored values.
  typedef T ValueType;

  /// Default constructor.
  Index() {}

  /// Constructor from const int array.
  Index( const T (&indexArray)[NDIM] )
  {
    for ( size_t i = 0; i < NDIM; ++i )
      this->m_Index[i] = indexArray[i];
  }

  /// Access operator.
  T& operator[]( const size_t idx )
  {
    return this->m_Index[idx];
  }

  /// Constant access operator.
  const T& operator[]( const size_t idx ) const
  {
    return this->m_Index[idx];
  }

  /// Get raw pointer to index array.
  operator T*()
  {
    return this->m_Index;
  }

  /// Get raw pointer-to-const to index array.
  operator const T*() const
  {
    return this->m_Index;
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
      this->m_Index[i] += rhs.m_Index[i];
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
      this->m_Index[i] -= rhs.m_Index[i];
    return *this;
  }
  
private:
  /// The actual index array.
  T m_Index[NDIM];
};

} // namespace cmtk

#endif // #ifndef __cmtkIndex_h_included_
