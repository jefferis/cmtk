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

#ifndef __cmtkMatrix4x4_h_included_
#define __cmtkMatrix4x4_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkFixedSquareMatrix.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkFixedVector.h>
#include <Base/cmtkMatrix3x3.h>

#include <System/cmtkConsole.h>
#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Homogeneous 4x4 transformation matrix.
template<class T=Types::Coordinate>
class Matrix4x4 :
    public FixedSquareMatrix<4,T>
{
public:
  /// This type instance.
  typedef Matrix4x4<T> Self;

  /// Base class..
  typedef FixedSquareMatrix<4,T> Superclass;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Identity transformation matrix.
  static const Self IdentityMatrix;

  /// Default constructor.
  Matrix4x4() {}

  /// Copy-from-baseclass constructor.
  Matrix4x4( const Superclass& other ) : Superclass( other ) {}

  /// Top left submatrix copy constructor.
  Matrix4x4( const Matrix3x3<T>& other );

  /** Array constructor.
   * If a NULL parameter is given, an uninitialized matrix is generated. This
   * is intended behaviour.
   */
  Matrix4x4( const T *const values ) : Superclass( values ) {}

  /// 2D array constructor.
  template<class T2> Matrix4x4( const T2 (&matrix)[4][4] ) : Superclass( matrix ) {}

  /// Compose from canonical parameters.
  Self& Compose( const Types::Coordinate params[15], const bool logScaleFactors = false );
  
  /// Decompose into affine parameters.
  bool Decompose( Types::Coordinate params[12], const Types::Coordinate *center = NULL, const bool logScaleFactors = false ) const;

  /** Change reference coordinate system.
   */
  Self& ChangeCoordinateSystem( const FixedVector<3,T>& newX, const  FixedVector<3,T>& newY );

  /// Return rotation around x axis.
  static Self RotateX( const T angle );
  
  /// Return rotation around y axis.
  static Self RotateY( const T angle );
  
  /// Return rotation around z axis.
  static Self RotateZ( const T angle );
};

template<typename T> const Matrix4x4<T> Matrix4x4<T>::IdentityMatrix;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMatrix4x4_h_included_
