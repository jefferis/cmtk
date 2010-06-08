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

#ifndef __cmtkVector3D_h_included_
#define __cmtkVector3D_h_included_

#include <cmtkconfig.h>

#include <cmtkFixedVector.h>

#include <math.h>

#include <cmtkMathUtil.h>
#include <cmtkTypes.h>
#include <cmtkConsole.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Vectors in 3D space.
 *@author Torsten Rohlfing
 *\deprecated This class will soon disappear and replaced with a simple typedef.
 */
class Vector3D 
  : public FixedVector<3,Types::Coordinate>
{
public:
  /// This class.
  typedef Vector3D Self;

  /// Parent class.
  typedef FixedVector<3,Types::Coordinate> Superclass;

  /// Create vector, but do not set any values.
  Vector3D () {}

  /// Create vector and initialize with a constant.
  explicit Vector3D ( const int v ) 
  {
    (*this)[0] = (*this)[1] = (*this)[2] = (Types::Coordinate) v;
  }

  /// Create vector from three single components
  explicit Vector3D ( const Types::Coordinate x, const Types::Coordinate y, const Types::Coordinate z ) 
  {
    (*this)[0] = x; (*this)[1] = y; (*this)[2] = z; 
  }

  /// Type conversion constructor template.
  template<class T2>
  Vector3D( const FixedVector<3,T2> rhs )
  {
    for ( size_t i = 0; i < 3; ++i )
      (*this)[i] = static_cast<Self::ValueType>( rhs[i] );
  }

  /// Set vector components.
  Vector3D& Set ( const Types::Coordinate x, const Types::Coordinate y, const Types::Coordinate z ) 
  {
    (*this)[0] = x; (*this)[1] = y; (*this)[2] = z;
    return *this;
  }

  /// Coordinatewise in place multiplication.
  static void CoordMultInPlace( Vector3D& p, const Vector3D& q );
  
  /// Coordinatewise divison.
  static const Vector3D CoordDiv( const Vector3D&, const Vector3D& );
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkVector3D_h_included_
