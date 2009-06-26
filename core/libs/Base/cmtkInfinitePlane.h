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

#ifndef __cmtkInfinitePlane_h_included_
#define __cmtkInfinitePlane_h_included_

#include <cmtkconfig.h>

#include <cmtkMacros.h>
#include <cmtkVector.h>
#include <cmtkVector3D.h>
#include <cmtkAffineXform.h>
#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for parameterized infinite planes.
 */
class InfinitePlane 
{
public:
  /// Default constructor.
  InfinitePlane();

  /** Origin of the coordinate system the plane parameters refer to.
   * This is NOT the origin of the plane. In fact, this coordinate does not
   * even need to be on the plane itself. It is merely the coordinate in space
   * relativ to which the plane is parameterized.
   */
  cmtkGetSetMacro(Vector3D,Origin);

  /// Set parameter Rho.
  void SetRho( const Types::Coordinate rho ) 
  {
    Rho = rho; this->Update();
  }
  
  /// Get parameter Rho.
  Types::Coordinate GetRho() const 
  { 
    return Rho; 
  }
  
  /// Set parameter Theta.
  void SetTheta( const Types::Coordinate theta ) 
  {
    Theta = theta; this->Update();
  }
  
  /// Get parameter Theta.
  Types::Coordinate GetTheta() const
  { 
    return Theta; 
  }
  
  /// Set parameter Phi.
  void SetPhi( const Types::Coordinate phi ) 
  {
    Phi = phi; this->Update();
  }
  
  /// Get parameter Phi.
  Types::Coordinate GetPhi() const 
  { 
    return Phi; 
  }

  /// Set normal vector and determine rotation angles from it.
  void SetNormal( const Vector3D& normal );

  /// Get all parameters.
  void GetParameters( CoordinateVector& v ) const 
  {
    v.SetDim( 6 );
    v[0] = Rho; v[1] = Theta; v[2] = Phi;
    v[3] = this->m_Origin[0]; v[4] = this->m_Origin[1]; v[5] = this->m_Origin[2];
  }
  
  /// Set all parameters.
  void SetParameters( const CoordinateVector& v ) 
  {
    Rho = v[0]; Theta = v[1]; Phi = v[2]; 
    this->m_Origin[0] = v[3]; this->m_Origin[1] = v[4]; this->m_Origin[2] = v[5];
    this->Update();
  }

  /** Determine which side of the plane a point is on.
   *@return 0, if given point is on the plane; +1 if point is on the same side
   * as Origin; -1, if point is on the other side, seen from Origin.
   */
  char GetWhichSide( const Vector3D& point ) const 
  {
    // move given origin to coordinate origin
    Vector3D p = point;
    p -= this->m_Origin;
    
    // compute line parameter of orthogonal projection of "point" onto this plane
    const Types::Coordinate intersect = Normal * p - Rho;
    return (intersect < 0) ? -1 : (intersect > 0) ? 1 : 0;
  }
  
  /// Mirror point with respect to plane.
  void Mirror( Vector3D& toPoint, const Vector3D& fromPoint ) const 
  {
    toPoint = fromPoint;
    this->MirrorInPlace( toPoint );
  }
  
  /// Mirror point in-place with respect to plane.
  void MirrorInPlace( Vector3D& point ) const 
  {
    // move given origin to coordinate origin
    point -= this->m_Origin;
    
    // compute line parameter of orthogonal projection of "point" onto
    // this plane and multiply by two to get parameter of mirrored point
    const Types::Coordinate intersect = 2 * (( Normal * point - Rho ) / SquareNormal);
    
    // compute mirrored point
    for ( int dim = 0; dim < 3; ++dim )
      point.XYZ[dim] -= intersect * Normal.XYZ[dim];

    // move given origin back to its given location
    point += this->m_Origin;
  }

  /// Project point onto plane.
  void Project( Vector3D& toPoint, const Vector3D& fromPoint ) const 
  {
    toPoint = fromPoint;
    this->ProjectInPlace( toPoint );
  }

  /// Project point onto plane in-place.
  void ProjectInPlace( Vector3D& point ) const 
  {
    // move given origin to coordinate origin
    point -= this->m_Origin;
    
    // compute line parameter of orthogonal projection of "point" onto
    // this plane
    const Types::Coordinate intersect = ( Normal * point - Rho ) / SquareNormal;
    
    // compute projected point
    for ( int dim = 0; dim < 3; ++dim )
      point.XYZ[dim] -= intersect * Normal.XYZ[dim];
    
    // move given origin back to its given location
    point += this->m_Origin;
  }

  /** Get transformation that aligns this plane with the coordinate system.
   * The object returned by this function represents a rigid transformation 
   * that aligns the normal vector of this object with one of the coordinate
   * axes.
   *@param axis The index of the axis to align the normal vector with (0, 1, or
   * 2 for x, y, or z, respectively).
   */
  AffineXform* GetAlignmentXform( const byte normalAxis = 0 ) const;

  /** Return normal vector.
   */
  const Vector3D& GetNormal() const { return Normal; }

private:
  /// Radius of tangent sphere with center at Origin.
  Types::Coordinate Rho;

  /// Rotation angle of tangent point around z axis through Origin.
  Types::Coordinate Theta;

  /// Elevation angle of tangent point over x/y plane through Origin.
  Types::Coordinate Phi;

  /// Plane normal.
  Vector3D Normal;

  /// Square norm of plane normal.
  Types::Coordinate SquareNormal;

  /// Update internal fields after change of parameters.
  void Update();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkInfinitePlane_h_included_
