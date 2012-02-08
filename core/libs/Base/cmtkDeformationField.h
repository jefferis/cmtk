/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#ifndef __cmtkDeformationField_h_included_
#define __cmtkDeformationField_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkWarpXform.h>
#include <Base/cmtkMacros.h>
#include <Base/cmtkVector.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkAffineXform.h>

#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Class for pixel-wise deformation field.
 *
 *\author Torsten Rohlfing
 */
class DeformationField : 
  /// Inherit generic grid-based nonrigid transformation interface.
  public WarpXform 
{
public:
  /// This class.
  typedef DeformationField Self;

  /// Parent class.
  typedef WarpXform Superclass;

  /// Smart pointer to DeformationField
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const DeformationField
  typedef SmartConstPointer<Self> SmartConstPtr;

public:
  /// Constructor.
  DeformationField( const UniformVolume* volume ) 
  {
    this->InitGrid( volume->Size, volume->m_Dims );
    this->m_Offset = volume->m_Offset;
  }
  
  /// Constructor.
  DeformationField( const FixedVector<3,Types::Coordinate>& domain, const DataGrid::IndexType& dims, const Types::Coordinate* offset = NULL ) 
  {
    this->InitGrid( domain, dims );
    if ( offset )
      {
      for ( int dim = 0; dim < 3; ++dim )
	this->m_Offset[dim] = offset[dim];
      }
  }
  
  /// Destructor.
  virtual ~DeformationField () {}

  /// Initialized internal data structures for new control point grid.
  virtual void InitGrid( const FixedVector<3,Types::Coordinate>& domain, const DataGrid::IndexType& dims )
  {
    this->Superclass::InitGrid( domain, dims );
    for ( int dim = 0; dim < 3; ++dim )
      {
      if ( dims[dim] > 1 )
	this->m_Spacing[dim] = domain[dim] / (dims[dim]-1);
      else
	this->m_Spacing[dim] = 1.0;
      this->m_InverseSpacing[dim] = 1.0 / this->m_Spacing[dim];
      }
    this->m_InverseAffineScaling[0] = this->m_InverseAffineScaling[1] = this->m_InverseAffineScaling[2] = this->m_GlobalScaling = 1.0;
  }
  
  /// Initialize control point positions, potentially with affine displacement.
  void InitControlPoints( const AffineXform* affineXform = NULL );

  /// Apply transformation to vector in-place.
  virtual void ApplyInPlace ( Self::SpaceVectorType& ) const;

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverse ( const Self::SpaceVectorType&, Self::SpaceVectorType&, const Types::Coordinate = 0.01  ) const 
  {
    // not implemented
    return false;
  }

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverseInPlace( Self::SpaceVectorType&, const Types::Coordinate = 0.01  ) const 
  {
    // not implemented
    return false;
  }

  /// Get a grid point from the deformed grid.
  virtual Self::SpaceVectorType GetTransformedGrid( const int idxX, const int idxY, const int idxZ ) const;
  
  /// Get a sequence of grid points from the deformed grid. 
  virtual void GetTransformedGridRow( Self::SpaceVectorType *const v, const int numPoints, const int idxX, const int idxY, const int idxZ ) const;
  
  /** Get the deformed position of a transformation control point.
   *\note This function does not return the shifted control point position,
   * but rather it applies the current transformation to the given control
   * point. It does so more efficiently than applying the transformation to
   * the explicit 3D coordinate of the control point, because most spline
   * coefficients vanish at control points and need not be considered.
   */
  virtual Self::SpaceVectorType GetDeformedControlPointPosition( const int, const int, const int ) const;
  
  /// Get Jacobian matrix.
  virtual void GetJacobian( const Self::SpaceVectorType& v, CoordinateMatrix3x3& J ) const;

  /// Compute Jacobian determinant at a certain location.
  virtual Types::Coordinate GetJacobianDeterminant ( const Self::SpaceVectorType& v ) const
  {
    CoordinateMatrix3x3 J;
    this->GetJacobian( v, J );
    return J.Determinant();
  }
  
  /// Return 1.0 since deformation field DOFs are always direct deformations in space units.
  virtual Types::Coordinate GetParamStep( const size_t, const Self::SpaceVectorType&, const Types::Coordinate mmStep = 1 ) const
  {
    return mmStep;
  }

protected:
  /** Clone transformation.
   *\todo This still needs to be implemented.
   */
  virtual Self* CloneVirtual () const { return NULL; }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeformationField_h_included_
