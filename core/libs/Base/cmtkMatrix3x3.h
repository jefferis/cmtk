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

#ifndef __cmtkMatrix3x3_h_included_
#define __cmtkMatrix3x3_h_included_

#include <cmtkconfig.h>

#include <cmtkTypes.h>
#include <cmtkFixedVector.h>
#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Homogeneous 3x3 transformation matrix.
template<class T=Types::Coordinate>
class Matrix3x3
{
public:
  /// This class.
  typedef Matrix3x3<T> Self;

  /** Default constructor: make identity matrix.
   *\note In order to create an uninitialized matrix (for speed), use
   * Matrix3x3<>( NULL ).
   *\todo We do not actually support shear yet.
   */
  Matrix3x3();

  /// Copy constructor.
  Matrix3x3( const Self& other );

  /** Array constructor.
   * If a NULL parameter is given, an uninitialized matrix is generated. This
   * is intended behaviour.
   */
  Matrix3x3( const T *const values ) 
  {
    if ( values ) this->Set( values );
  }

  /// 2D array constructor.
  template<class T2>
  Matrix3x3( const T2 (&matrix)[3][3] ) 
  {
    for ( size_t j = 0; j < 3; ++j )
      for ( size_t i = 0; i < 3; ++i )
	this->Matrix[j][i] = matrix[j][i];
  }

  /// Set from array of entries.
  Self& Set( const T *const values );

  /// Set to constant value.
  Self& Fill( const T value );

  /// Inversion operator (in place).
  Self& Invert2x2();

  /// Inversion operator (in place) as a 3D non-homogeneous matrix.
  Self& Invert3x3();

  /// Transpose operator.
  Self GetTranspose() const;

  /// Index operator.
  T* operator[]( const size_t i )
  { 
    return this->Matrix[i]; 
  }

  /// Constant index operator.
  const T* operator[]( const size_t i ) const
  { 
    return this->Matrix[i]; 
  }
  
  /// Compose from canonical parameters.
  Self& Compose( const Types::Coordinate params[8] );
  
  /// Decompose into affine parameters.
  bool Decompose( Types::Coordinate params[8], const Types::Coordinate *center = NULL ) const;

  /// In-place multiplication operator.
  Self& operator*=( const Self& other );
  
  /// In-place scalar multiplication operator.
  Self& operator*=( const T scalar );
  
  /// Multiplication operator.
  const Self operator*( const Self& other ) const;

  /// Assignment operator.
  Self& operator=( const Self& other );

  /* Multiply with homogeneous 3d vector (implicitly made homogeneous).
   *\note The '&' declaration of both arguments forces the C++ compiler to
   * maintain them as explicitly sized arrays, rather than collapsing to
   * pointers. This is what allows us to overload this function based on the
   * array size (!) of its arguments.
   */
  template<class T2> void Multiply( const FixedVector<2,T2>& u, FixedVector<2,T2>& v ) const 
  {
    for ( int idx=0; idx < 2; ++idx ) 
      v[idx] = u[0]*Matrix[0][idx] + u[1]*Matrix[1][idx] + Matrix[2][idx];
  }
  
  /** Multiply with actual 3d vector.
   *\note The '&' declaration of both arguments forces the C++ compiler to
   * maintain them as explicitly sized arrays, rather than collapsing to
   * pointers. This is what allows us to overload this function based on the
   * array size (!) of its arguments.
   */
  template<class T2> void Multiply( const FixedVector<3,T2>& u,  FixedVector<3,T2>& v ) const 
  {
    for ( int idx=0; idx < 3; ++idx ) 
      v[idx] = u[0]*Matrix[0][idx] + u[1]*Matrix[1][idx] + u[2]*Matrix[2][idx];
  }
  
  /// Multiply in place with 3d vector (will implicitly be made homogeneous).
  template<class T2> void Multiply( FixedVector<2,T2>& v ) const 
  {
    const FixedVector<2,T2> u( v );
    this->Multiply( u, v );
  }
  
  /// Multiply in place with actual 3d vector.
  template<class T2> void Multiply( FixedVector<3,T2>& v ) const 
  {
    const FixedVector<3,T2> u( v );
    this->Multiply( u, v );
  }

  /// Get determinant.
  T Determinant() const 
  {
    return ( this->Matrix[0][0]*this->Matrix[1][1]*this->Matrix[2][2] + 
	     this->Matrix[0][1]*this->Matrix[1][2]*this->Matrix[2][0] + 
	     this->Matrix[0][2]*this->Matrix[1][0]*this->Matrix[2][1] - 
	     this->Matrix[0][2]*this->Matrix[1][1]*this->Matrix[2][0] - 
	     this->Matrix[0][0]*this->Matrix[1][2]*this->Matrix[2][1] - 
	     this->Matrix[0][1]*this->Matrix[1][0]*this->Matrix[2][2] );
  }

  /// Get Frobenius norm.
  T FrobeniusNorm() const;

  /// Compute eigenvalues.
  void ComputeEigenvalues( T (&lambda)[3] ) const;

private:
  /// The actual matrix.
  T Matrix[3][3];
};

/// Define coordinate matrix.
typedef Matrix3x3<Types::Coordinate> CoordinateMatrix3x3;

/// Output object to console.
template<class T>
inline
Console& operator<< ( Console& stream, const Matrix3x3<T>& m )
{
  stream << "3x3 Matrix:\n";
  for ( int i = 0; i < 3; ++i ) 
    {
    for ( int j = 0; j < 3; ++j )
      stream << m[i][j] << "\t";
    stream << "\n";
    }
  return stream;
}

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMatrix3x3_h_included_
