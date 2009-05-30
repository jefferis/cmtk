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

#include <cmtkVector3D.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

Vector3D& 
Vector3D::SetNormal( const Vector3D& x, const Vector3D& y )
{
  this->Set( x.XYZ[1] * y.XYZ[2] - x.XYZ[2] * y.XYZ[1], x.XYZ[2] * y.XYZ[0] - x.XYZ[0] * y.XYZ[2], x.XYZ[0] * y.XYZ[1] - x.XYZ[1] * y.XYZ[0] );

  return *this;
}

byte
Vector3D::GetMaxComponentDimension() const
{
  if ( (XYZ[0] > XYZ[1]) && (XYZ[0] > XYZ[2]) )
    return AXIS_X;

  if ( (XYZ[1] > XYZ[2]) )
    return AXIS_Y;

  return AXIS_Z;
}

Vector3D
operator+ ( const Vector3D& p, const Vector3D& delta )
{
  return Vector3D( p[0]+delta[0], p[1]+delta[1], p[2]+delta[2] );
}

Vector3D
operator- ( const Vector3D& p, const Vector3D& delta )
{
  return Vector3D( p[0]-delta[0], p[1]-delta[1], p[2]-delta[2] );
}

Vector3D
operator* ( const Types::Coordinate c, const Vector3D& p ) 
{
  return Vector3D( c*p[0], c*p[1], c*p[2] );
}

Vector3D 
Vector3D::CoordMult ( const Vector3D& p, const Vector3D& q ) 
{
  return Vector3D( p[0]*q[0], p[1]*q[1], p[2]*q[2]);
}

void
Vector3D::CoordMultInPlace( Vector3D& p, const Vector3D& q ) 
{
  for ( int dim = 0; dim < 3; ++dim ) p.XYZ[dim] *= q.XYZ[dim];
}

Vector3D 
Vector3D::CoordDiv ( const Vector3D& p, const Vector3D& q ) 
{
  return Vector3D( p[0]/q[0], p[1]/q[1], p[2]/q[2]);
}

Types::Coordinate
operator* ( const Vector3D& p, const Vector3D& q ) 
{
  return p[0]*q[0]+p[1]*q[1]+p[2]*q[2];
}

} // namespace cmtk
