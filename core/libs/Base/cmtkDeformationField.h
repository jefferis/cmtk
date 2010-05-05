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

#ifndef __cmtkDeformationField_h_included_
#define __cmtkDeformationField_h_included_

#include <cmtkconfig.h>

#include <cmtkWarpXform.h>

#include <cmtkMacros.h>
#include <cmtkVector.h>
#include <cmtkVector3D.h>

#include <cmtkUniformVolume.h>

#include <cmtkAffineXform.h>
#include <cmtkSmartPtr.h>

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

public:
  /// Constructor.
  DeformationField( const UniformVolume* volume ) 
  {
    this->InitGrid( volume->Size, volume->m_Dims );
    this->m_Offset = volume->m_Offset;
  }
  
  /// Constructor.
  DeformationField( const Types::Coordinate domain[3], const DataGrid::IndexType& dims, const Types::Coordinate* offset = NULL ) 
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

  /** Clone transformation.
   *\todo This still needs to be implemented.
   */
  virtual Self* Clone () const { return NULL; }

  /// Initialized internal data structures for new control point grid.
  virtual void InitGrid( const Types::Coordinate domain[3], const DataGrid::IndexType& dims )
  {
    this->Superclass::InitGrid( domain, dims );
    for ( int dim = 0; dim < 3; ++dim )
      {
      if ( dims[dim] > 1 )
	this->Spacing[dim] = domain[dim] / (dims[dim]-1);
      else
	this->Spacing[dim] = 1.0;
      this->InverseSpacing[dim] = 1.0 / this->Spacing[dim];
      }
    this->InverseAffineScaling[0] = this->InverseAffineScaling[1] = this->InverseAffineScaling[2] = this->GlobalScaling = 1.0;
  }
  
  /// Initialize control point positions, potentially with affine displacement.
  void InitControlPoints( const AffineXform* affineXform = NULL );

  /// Apply transformation to vector in-place.
  virtual void ApplyInPlace ( Vector3D& ) const;

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverse ( const Vector3D&, Vector3D&, const Types::Coordinate = 0.01  ) const 
  {
    // not implemented
    return false;
  }

  /** Return origin of warped vector.
   */
  virtual bool ApplyInverseInPlace( Vector3D&, const Types::Coordinate = 0.01  ) const 
  {
    // not implemented
    return false;
  }

  /// Get a grid point from the deformed grid.
  virtual void GetTransformedGrid( Vector3D& v, const int idxX, const int idxY, const int idxZ ) const;
  
  /// Get a sequence of grid points from the deformed grid. 
  virtual void GetTransformedGridSequence( Vector3D *const v, const int numPoints, const int idxX, const int idxY, const int idxZ ) const;
  
  /// Get Jacobian matrix.
  virtual void GetJacobian( const Vector3D& v, CoordinateMatrix3x3& J ) const;

  /// Compute Jacobian determinant at a certain location.
  virtual Types::Coordinate GetJacobianDeterminant ( const Vector3D& v ) const
  {
    CoordinateMatrix3x3 J;
    this->GetJacobian( v, J );
    return J.Determinant();
  }
  
  /// Return 1.0 since deformation field DOFs are always direct deformations in space units.
  virtual Types::Coordinate GetParamStep( const size_t, const Types::Coordinate*, const Types::Coordinate = 1 ) const
  {
    return 1.0;
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeformationField_h_included_
