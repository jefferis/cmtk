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

  /// Create vector from array of double precision numbers
  explicit Vector3D ( const double *array ) 
  { 
    (*this)[0]=(Types::Coordinate)array[0]; 
    (*this)[1]=(Types::Coordinate)array[1]; 
    (*this)[2]=(Types::Coordinate)array[2]; 
  }

  /// Create vector from array of floats
  explicit Vector3D ( const float *array ) 
  { 
    (*this)[0]=(Types::Coordinate)array[0]; 
    (*this)[1]=(Types::Coordinate)array[1]; 
    (*this)[2]=(Types::Coordinate)array[2]; 
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

  /// Set this vector to be the normal vector of two other vectors.
  Vector3D& SetNormal( const Vector3D& x, const Vector3D& y );

  /// Calculate maximum vector norm.
  Types::Coordinate MaxNorm () const 
  { 
    return (Types::Coordinate) std::max( fabs((*this)[0]), std::max( fabs((*this)[1]), fabs((*this)[2]) ) );
  }
  
  /// Get dimension with maximum vector component.
  byte GetMaxComponentDimension() const;

  /// Test if vector is all zero.
  bool IsNull () const 
  { 
    for ( int dim=0; dim<3; ++dim )
      if ( (*this)[dim] != 0 ) return false;
    return true;
  }
  
  /// Coordinatewise multiplication.
  static const Vector3D CoordMult( const Vector3D&, const Vector3D& );

  /// Coordinatewise in place multiplication.
  static void CoordMultInPlace( Vector3D& p, const Vector3D& q );
  
  /// Coordinatewise divison.
  static const Vector3D CoordDiv( const Vector3D&, const Vector3D& );
};

inline Types::Coordinate
operator* ( const Vector3D& p, const Vector3D& q ) 
{
  return p[0]*q[0]+p[1]*q[1]+p[2]*q[2];
}

/// Stream output operator.
inline
Console& operator<< ( Console& stream, const Vector3D& v )
{
  for ( int i = 0; i < 3; ++i )
    stream << v[i] << "\t";
  stream << "\n";
  return stream;
}

//@}

} // namespace cmtk

#endif // #ifdef __cmtkVector3D_h_included_
