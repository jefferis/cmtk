/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkParametricPlane_h_included_
#define __cmtkParametricPlane_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMacros.h>
#include <Base/cmtkVector.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkMathUtil.h>
#include <Base/cmtkUnits.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for parameterized infinite planes.
 */
class ParametricPlane 
{
public:
  /// This class.
  typedef ParametricPlane Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Coordinate vector.
  typedef FixedVector<3,Types::Coordinate> CoordinateVectorType;

  /// Default constructor.
  ParametricPlane();

  /** Origin of the coordinate system the plane parameters refer to.
   * This is NOT the origin of the plane. In fact, this coordinate does not
   * even need to be on the plane itself. It is merely the coordinate in space
   * relative to which the plane is parameterized.
   */
  cmtkGetSetMacro(Self::CoordinateVectorType,Origin);

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
  void SetTheta( const Units::Degrees& theta ) 
  {
    Theta = theta; this->Update();
  }
  
  /// Get parameter Theta.
  const Units::Degrees GetTheta() const
  { 
    return Theta; 
  }
  
  /// Set parameter Phi.
  void SetPhi( const Units::Degrees& phi ) 
  {
    Phi = phi; this->Update();
  }
  
  /// Get parameter Phi.
  const Units::Degrees GetPhi() const 
  { 
    return Phi; 
  }

  /// Set normal vector and determine rotation angles from it.
  void SetNormal( const Self::CoordinateVectorType& normal );

  /// Get all parameters.
  void GetParameters( CoordinateVector& v ) const 
  {
    v.SetDim( 6 );
    v[0] = Rho; v[1] = Theta.Value(); v[2] = Phi.Value();
    v[3] = this->m_Origin[0]; v[4] = this->m_Origin[1]; v[5] = this->m_Origin[2];
  }
  
  /// Set all parameters.
  void SetParameters( const CoordinateVector& v ) 
  {
    Rho = v[0]; Theta = Units::Degrees( v[1] ); Phi = Units::Degrees( v[2] ); 
    this->m_Origin[0] = v[3]; this->m_Origin[1] = v[4]; this->m_Origin[2] = v[5];
    this->Update();
  }

  /** Determine which side of the plane a point is on.
   *@return 0, if given point is on the plane; +1 if point is on the same side
   * as Origin; -1, if point is on the other side, seen from Origin.
   */
  char GetWhichSide( const Self::CoordinateVectorType& point ) const 
  {
    // move given origin to coordinate origin
    Self::CoordinateVectorType p = point;
    p -= this->m_Origin;
    
    // compute line parameter of orthogonal projection of "point" onto this plane
    const Types::Coordinate intersect = Normal*p - Rho;
    return (intersect < 0) ? -1 : (intersect > 0) ? 1 : 0;
  }
  
  /// Mirror point with respect to plane.
  void Mirror( Self::CoordinateVectorType& toPoint, const Self::CoordinateVectorType& fromPoint ) const 
  {
    toPoint = fromPoint;
    this->MirrorInPlace( toPoint );
  }
  
  /// Mirror point in-place with respect to plane.
  void MirrorInPlace( Self::CoordinateVectorType& point ) const 
  {
    // move given origin to coordinate origin
    point -= this->m_Origin;
    
    // compute line parameter of orthogonal projection of "point" onto
    // this plane and multiply by two to get parameter of mirrored point
    const Types::Coordinate intersect = 2 * (( Normal * point - Rho ) / SquareNormal);
    
    // compute mirrored point
    for ( int dim = 0; dim < 3; ++dim )
      point[dim] -= intersect * Normal[dim];

    // move given origin back to its given location
    point += this->m_Origin;
  }

  /// Project point onto plane.
  void Project( Self::CoordinateVectorType& toPoint, const Self::CoordinateVectorType& fromPoint ) const 
  {
    toPoint = fromPoint;
    this->ProjectInPlace( toPoint );
  }

  /// Project point onto plane in-place.
  void ProjectInPlace( Self::CoordinateVectorType& point ) const 
  {
    // move given origin to coordinate origin
    point -= this->m_Origin;
    
    // compute line parameter of orthogonal projection of "point" onto
    // this plane
    const Types::Coordinate intersect = ( Normal * point - Rho ) / SquareNormal;
    
    // compute projected point
    for ( int dim = 0; dim < 3; ++dim )
      point[dim] -= intersect * Normal[dim];
    
    // move given origin back to its given location
    point += this->m_Origin;
  }

  /** Get transformation that aligns this plane with the coordinate system.
   * The object returned by this function represents a rigid transformation 
   * that aligns the normal vector of this object with one of the coordinate
   * axes.
   *@param normalAxis The index of the axis to align the normal vector with (0, 1, or
   * 2 for x, y, or z, respectively).
   */
  AffineXform* GetAlignmentXform( const byte normalAxis = 0 ) const;

  /// Get affine transformation matrix of the mirror transform.
  AffineXform::MatrixType GetMirrorXformMatrix() const;

  /** Return normal vector.
   */
  const Self::CoordinateVectorType& GetNormal() const { return Normal; }

private:
  /// Radius of tangent sphere with center at Origin.
  Types::Coordinate Rho;

  /// Rotation angle of tangent point around z axis through Origin.
  Units::Degrees Theta;

  /// Elevation angle of tangent point over x/y plane through Origin.
  Units::Degrees Phi;

  /// Plane normal.
  Self::CoordinateVectorType Normal;

  /// Square norm of plane normal.
  Types::Coordinate SquareNormal;

  /// Update internal fields after change of parameters.
  void Update();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkParametricPlane_h_included_
