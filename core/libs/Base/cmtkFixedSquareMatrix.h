/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkFixedSquareMatrix_h_included_
#define __cmtkFixedSquareMatrix_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkFixedVector.h>

#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Fixed-size rank-2 matrix.
template<size_t NDIM,class TSCALAR=Types::Coordinate>
class FixedSquareMatrix
{
public:
  /// This class.
  typedef FixedSquareMatrix<NDIM,TSCALAR> Self;

  /// The scalar data type.
  typedef TSCALAR ScalarType;

  /// The matrix dimension.
  static const size_t Dimension = NDIM;

  /** Default constructor: make identity matrix.
   *\note In order to create an uninitialized matrix (for speed), use
   * FixedSquareMatrix<>( NULL ).
   *\todo We do not actually support shear yet.
   */
  FixedSquareMatrix();

  /// Copy constructor.
  FixedSquareMatrix( const Self& other );

  /** Array constructor.
   * If a NULL parameter is given, an uninitialized matrix is generated. This
   * is intended behaviour.
   */
  FixedSquareMatrix( const typename Self::ScalarType *const values ) 
  {
    if ( values ) 
      this->Set( values );
  }

  /// 2D array constructor.
  template<class T2>
  FixedSquareMatrix( const T2 (&matrix)[NDIM][NDIM] ) 
  {
    for ( size_t j = 0; j < NDIM; ++j )
      for ( size_t i = 0; i < NDIM; ++i )
	this->m_Matrix[j][i] = matrix[j][i];
  }

  /// Set from array of entries.
  Self& Set( const typename Self::ScalarType *const values );

  /// Set to constant value.
  Self& Fill( const typename Self::ScalarType& value );

  /// Get inverse matrix.
  Self GetInverse() const;

  /// Get transpose matrix.
  Self GetTranspose() const;

  /// Index operator.
  typename Self::ScalarType* operator[]( const size_t i )
  { 
    return this->m_Matrix[i]; 
  }

  /// Constant index operator.
  const typename Self::ScalarType* operator[]( const size_t i ) const
  { 
    return this->m_Matrix[i]; 
  }
  
  /// In-place multiplication operator.
  Self& operator*=( const Self& other );
  
  /// In-place scalar multiplication operator.
  Self& operator*=( const typename Self::ScalarType& scalar );
  
  /// Multiplication operator.
  const Self operator*( const Self& other ) const;

  /// Assignment operator.
  Self& operator=( const Self& other );

  /// Get Frobenius norm.
  typename Self::ScalarType FrobeniusNorm() const;

protected:
  /// The actual matrix.
  typename Self::ScalarType m_Matrix[NDIM][NDIM];
};

/// Multiplication with vector operator.
template<size_t NDIM,class T> 
FixedVector<NDIM,T>
operator*( const FixedVector<NDIM,T>& u, const FixedSquareMatrix<NDIM,T>& M )
{
  FixedVector<NDIM,T> v;
  for ( size_t i = 0; i<NDIM; ++i ) 
    {
    v[i] = u[0]*M[0][i];
    for ( size_t j = 1; j<NDIM; ++j ) 
      {
      v[i] += u[j]*M[j][i];
      }
    }
  return v;
}

/// In-place vector-matrix multiplication operation.
template<size_t NDIM,class T> 
FixedVector<NDIM,T>&
operator*=( FixedVector<NDIM,T>& u, const FixedSquareMatrix<NDIM,T>& M )
{
  return u = u*M;
}

/// Multiplication with homogeneous vector operator.
template<size_t NDIM,class T> 
FixedVector<NDIM,T>
operator*( const FixedVector<NDIM,T>& u, const FixedSquareMatrix<NDIM+1,T>& M )
{
  FixedVector<NDIM,T> v;
  for ( size_t i = 0; i<NDIM; ++i ) 
    {
    v[i] = u[0]*M[0][i];
    for ( size_t j = 1; j<NDIM; ++j ) 
      {
      v[i] += u[j]*M[j][i];
      }
    // add 1x last column element for implicitly homogeneous vector.
    v[i] += M[NDIM][i];
    }
  return v;
}

/// In-place homogeneous multiplication with vector operation.
template<size_t NDIM,class T> 
FixedVector<NDIM,T>&
operator*=( FixedVector<NDIM,T>& u, const FixedSquareMatrix<NDIM+1,T>& M )
{
  return u = u*M;
}

/// Output object to console.
template<size_t NDIM,class TSCALAR>
inline
Console& operator<< ( Console& stream, const FixedSquareMatrix<NDIM,TSCALAR>& m )
{
  stream << NDIM << "x" << NDIM << " Matrix:\n";
  for ( int i = 0; i < NDIM; ++i ) 
    {
    for ( int j = 0; j < NDIM; ++j )
      stream << m[i][j] << "\t";
    stream << "\n";
    }
  return stream;
}

//@}

} // namespace cmtk

#include "cmtkFixedSquareMatrix.txx"

#endif // #ifndef __cmtkFixedSquareMatrix_h_included_
