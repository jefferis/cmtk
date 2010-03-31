/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <stdio.h>
#include <string.h>
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
 */
class Vector3D 
{
public:
  /// This class.
  typedef Vector3D Self;

  /**@name Vector elements */
  Types::Coordinate XYZ[3];

  /// Create vector, but do not set any values.
  Vector3D () {}

  /// Create vector and initialize with a constant.
  Vector3D ( const int v ) 
  {
    XYZ[0] = XYZ[1] = XYZ[2] = (Types::Coordinate) v;
  }

  /// Create vector from array of double precision numbers
  Vector3D ( const double *array ) 
  { 
    XYZ[0]=(Types::Coordinate)array[0]; 
    XYZ[1]=(Types::Coordinate)array[1]; 
    XYZ[2]=(Types::Coordinate)array[2]; 
  }

  /// Create vector from array of floats
  Vector3D ( const float *array ) 
  { 
    XYZ[0]=(Types::Coordinate)array[0]; 
    XYZ[1]=(Types::Coordinate)array[1]; 
    XYZ[2]=(Types::Coordinate)array[2]; 
  }

  /// Create vector from three single components
  Vector3D ( const Types::Coordinate x, const Types::Coordinate y, const Types::Coordinate z ) 
  {
    XYZ[0] = x; XYZ[1] = y; XYZ[2] = z; 
  }

  /// Clone vector
  Vector3D ( const Vector3D& other ) 
  {
    XYZ[0] = other.XYZ[0]; XYZ[1] = other.XYZ[1]; XYZ[2] = other.XYZ[2];
  }

  /// Set x-component of this vector.
  void SetX ( const Types::Coordinate value ) { XYZ[0]=value; }
 
  /// Set y-component of this vector.
  void SetY ( const Types::Coordinate value ) { XYZ[1]=value; }

  /// Set z-component of this vector.
  void SetZ ( const Types::Coordinate value ) { XYZ[2]=value; }

  /// Return x-component of this vector.
  Types::Coordinate X() const { return XYZ[0]; }

  /// Return y-component of this vector.
  Types::Coordinate Y() const { return XYZ[1]; }

  /// Return z-component of this vector.
  Types::Coordinate Z() const { return XYZ[2]; }

  /// Set vector components.
  Vector3D& Set ( const Types::Coordinate x, const Types::Coordinate y, const Types::Coordinate z ) 
  {
    XYZ[0] = x; XYZ[1] = y; XYZ[2] = z;
    return *this;
  }

  /// Set vector components.
  Vector3D& Set ( const double* array ) 
  {
    for ( int dim=0; dim<3; ++dim )
      XYZ[dim] = (Types::Coordinate) array[dim];
    return *this;
  }

  /// Set vector components.
  Vector3D& Set ( const float* array ) 
  {
    for ( int dim=0; dim<3; ++dim ) 
      XYZ[dim] = (Types::Coordinate) array[dim];
    return *this;
  }

  /// Set single coordinate by index.
  Vector3D& Set ( const int index, const Types::Coordinate value ) {
    XYZ[index] = value;
    return *this;
  }

  /// Set this vector to be the normal vector of two other vectors.
  Vector3D& SetNormal( const Vector3D& x, const Vector3D& y );

  /// Vector assignment.
  Vector3D& operator = ( const Vector3D& other ) 
  {
    memcpy( XYZ, other.XYZ, sizeof( XYZ ) ); 
    return *this; 
  }

  /// Test for vector equality.
  bool operator== ( const Vector3D& other ) const 
  { 
    for ( int dim=0; dim<3; ++dim )
      if ( XYZ[dim] != other[dim] ) return false;
    return true;
  }
  
  /// Test for vector inequality.
  bool operator!= ( const Vector3D& other ) const 
  {
    for ( int dim=0; dim<3; ++dim )
      if ( XYZ[dim] != other[dim] ) return true;
    return false;
  }
  
  /// Test for "smaller or equal".
  bool operator<= ( const Vector3D& other ) const 
  {
    for ( int dim=0; dim<3; ++dim )
      if ( XYZ[dim] > other[dim] ) return false;
    return true;
  }

  /// Test for "greater or equal".
  bool operator>= ( const Vector3D& other ) const 
  {
    for ( int dim=0; dim<3; ++dim )
      if ( XYZ[dim] < other[dim] ) return false;
    return true;
  }

  /// Calculate maximum vector norm.
  Types::Coordinate MaxNorm () const 
  { 
    return (Types::Coordinate) std::max( fabs(XYZ[0]), std::max( fabs(XYZ[1]), fabs(XYZ[2]) ) );
  }
  
  /// Calculate square of vector (that is, square of Euclid's norm).
  Types::Coordinate Square () const 
  { 
    return (Types::Coordinate)( MathUtil::Square( XYZ[0] ) + MathUtil::Square( XYZ[1] ) + MathUtil::Square( XYZ[2] ) ); 
  }

  /// Calculate square of Euclid distance.
  static Types::Coordinate SquareEuclidDistance( const Self& x, const Self& y )
  {
    return MathUtil::Square(x.XYZ[0]-y.XYZ[0]) + MathUtil::Square(x.XYZ[1]-y.XYZ[1]) + MathUtil::Square(x.XYZ[2]-y.XYZ[2]);
  }

  /// Calculate Euclid's vector norm.
  Types::Coordinate EuclidNorm () const 
  { 
    return static_cast<Types::Coordinate>( sqrt( this->Square() ) ); 
  }
  
  /// Get dimension with maximum vector component.
  byte GetMaxComponentDimension() const;

  /// Test if vector is all zero.
  bool IsNull () const 
  { 
    for ( int dim=0; dim<3; ++dim )
      if ( XYZ[dim] != 0 ) return false;
    return true;
  }
  
  /// Get reference to constant vector element by coordinate index.
  inline const Types::Coordinate& operator[] ( const int index ) const 
  {
    return this->XYZ[index];
  }
  
  /// Get reference to non-constant vector element by coordinate index.
  inline Types::Coordinate& operator[] ( const int index ) 
  {
    return this->XYZ[index];
  }

  /// Increment vector by another.
  Vector3D& operator+= ( const Vector3D& other )
  {
    for ( int dim=0; dim<3; ++dim ) 
      XYZ[dim] += other[dim];
    return *this;
  }

  /// Decrement vector by another.
  Vector3D& operator-= ( const Vector3D& other ) 
  {
    for ( int dim=0; dim<3; ++dim ) 
      XYZ[dim] -= other[dim];
    return *this;
  }

  /// Multiply by a scalar.
  Vector3D& operator*= ( const Types::Coordinate a ) 
  {
    for ( int dim=0; dim<3; ++dim ) 
      XYZ[dim] *= a;
    return *this;
  }

  /// Divide by a scalar.
  Vector3D& operator/= ( const Types::Coordinate a ) 
  {
    Types::Coordinate inv_a = (Types::Coordinate) (1.0/a);
    for ( int dim=0; dim<3; ++dim ) 
      XYZ[dim] *= inv_a;
    return *this;
  }

  /// Unary minus.
  Vector3D operator-() { return Vector3D(-XYZ[0],-XYZ[1],-XYZ[2]); }

  /**@name Friend operators */
  //@{
  /// Add two vectors
  friend const Vector3D operator + ( const Vector3D&, const Vector3D& );

  /// Subtract one vector from another
  friend const Vector3D operator - ( const Vector3D&, const Vector3D& );

  /// Inner product
  friend Types::Coordinate operator * ( const Vector3D&, const Vector3D& );

  /// Scalar multiplication
  friend const Vector3D operator * ( const Types::Coordinate, const Vector3D& );

  /// Coordinatewise multiplication.
  static const Vector3D CoordMult( const Vector3D&, const Vector3D& );

  /// Coordinatewise in place multiplication.
  static void CoordMultInPlace( Vector3D& p, const Vector3D& q );

  /// Coordinatewise divison.
  static const Vector3D CoordDiv( const Vector3D&, const Vector3D& );
  //@}
};

const Vector3D operator * ( const Types::Coordinate, const Vector3D& );

/// Stream output operator.
inline
Console& operator<< ( Console& stream, const Vector3D& v )
{
  for ( int i = 0; i < 3; ++i )
    stream << v.XYZ[i] << "\t";
  stream << "\n";
  return stream;
}

//@}

} // namespace cmtk

#endif // #ifdef __cmtkVector3D_h_included_
