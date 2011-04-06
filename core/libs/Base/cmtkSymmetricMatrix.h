/*
//
//  Copyright 2011 SRI International
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

#ifndef __cmtkSymmetricMatrix_h_included_
#define __cmtkSymmetricMatrix_h_included_

#include <cmtkconfig.h>

#include <vector>
#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class template for variable-size symmetric matrix.
template<class TElement>
class SymmetricMatrix
{
public:
  /// This class.
  typedef SymmetricMatrix<TElement> Self;

  /// Constructor.
  SymmetricMatrix( const size_t dim = 0 ) : m_Dim( dim ), m_MatrixElements( Self::NumberOfElements( dim ) ) {}

  /// Reference to matrix element.
  TElement& operator()( const size_t i, const size_t j )
  {
    return this->m_MatrixElements[std::max(i,j) + std::min(i,j)*this->m_Dim];
  }
  
  /// Reference to const matrix element.
  const TElement& operator()( const size_t i, const size_t j ) const
  {
    return this->m_MatrixElements[std::max(i,j) + std::min(i,j)*this->m_Dim];
  }

  /// Get matrix dimension.
  size_t Dim() const
  {
    return this->m_Dim;
  }
  
  /// Resize matrix.
  void Resize( const size_t newDim )
  {
    this->m_Dim = newDim;
    this->m_MatrixElements.resize( Self::NumberOfElements( newDim ) );
  }
  
  /// Resize matrix with explicit initializer for newly allocated elements.
  void Resize( const size_t newDim, const TElement initValue )
  {
    this->m_Dim = newDim;
    this->m_MatrixElements.resize( Self::NumberOfElements( newDim ), initValue );
  }

  /// Equality operator.
  bool operator==( const Self& other ) const
  {
    if ( this->m_Dim != other.m_Dim )
      return false;
    return (this->m_MatrixElements == other.m_MatrixElements);
  }
  
  /// Inequality operator.
  bool operator!=( const Self& other ) const
  {
    return ! (*this == other);
  }
  
private:
  /// Static member: compute number of elements from dimension.
  static size_t NumberOfElements( const size_t dim )
  {
    return dim*(dim+1) / 2;
  }

  /// Matrix dimension.
  size_t m_Dim;

  /// Vector that stores matrix elements.
  std::vector<TElement> m_MatrixElements;
};

//@}

} // namespace cmtk


#endif // #ifndef __cmtkSymmetricMatrix_h_included_
