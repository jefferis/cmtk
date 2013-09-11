/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include <Base/cmtkFixedSquareMatrix.h>
#include <Base/cmtkFixedVector.h>

#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Homogeneous 3x3 transformation matrix.
template<class T=Types::Coordinate>
class Matrix3x3 :
    public FixedSquareMatrix<3,T>
{
public:
  /// This class.
  typedef Matrix3x3<T> Self;

  /// The floating point element type.
  typedef T ElementType;

  /// Parent class.
  typedef FixedSquareMatrix<3,T> Superclass;

  /// Default constructor.
  Matrix3x3() {}

  /// Copy constructor.
  Matrix3x3( const Superclass& other ) : Superclass( other ) {}

  /** Array constructor.
   * If a NULL parameter is given, an uninitialized matrix is generated. This
   * is intended behaviour.
   */
  Matrix3x3( const typename Self::ElementType *const values ) : Superclass( values ) {}
  
  /// 2D array constructor.
  template<class T2> Matrix3x3( const T2 (&matrix)[3][3] ) : Superclass( matrix ) {}

  /// Compose from canonical parameters.
  Self& Compose( const typename Self::ElementType params[8] );
  
  /// Decompose into affine parameters.
  bool Decompose( Self::ElementType params[8], const typename Self::ElementType *center = NULL ) const;

  /// Get determinant.
  typename Self::ElementType Determinant() const 
  {
    return ( (*this)[0][0]*(*this)[1][1]*(*this)[2][2] + 
	     (*this)[0][1]*(*this)[1][2]*(*this)[2][0] + 
	     (*this)[0][2]*(*this)[1][0]*(*this)[2][1] - 
	     (*this)[0][2]*(*this)[1][1]*(*this)[2][0] - 
	     (*this)[0][0]*(*this)[1][2]*(*this)[2][1] - 
	     (*this)[0][1]*(*this)[1][0]*(*this)[2][2] );
  }

  /// Compute eigenvalues.
  void ComputeEigenvalues( typename Self::ElementType (&lambda)[3] ) const;
};

/// Define coordinate matrix.
typedef Matrix3x3<Types::Coordinate> CoordinateMatrix3x3;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMatrix3x3_h_included_
