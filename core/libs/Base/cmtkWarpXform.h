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

#ifndef __cmtkWarpXform_h_included_
#define __cmtkWarpXform_h_included_

#include <cmtkconfig.h>

#include <cmtkXform.h>

#include <cmtkMacros.h>
#include <cmtkVector.h>
#include <cmtkVector3D.h>
#include <cmtkRect3D.h>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>

#include <cmtkAffineXform.h>
#include <cmtkBitVector.h>

#include <cmtkControlPoint.h>
#include <cmtkMatchedLandmarkList.h>

#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Common base class for free-form-deformation-based warps.
 */
class WarpXform : 
  /// Inherit generic transformation interface.
  public Xform 
{
public:
  /// Smart pointer to WarpXform
  typedef SmartPointer<WarpXform> SmartPtr;

  /// Smart pointer to const WarpXform
  typedef SmartConstPointer<WarpXform> SmartConstPtr;

  /// Dimensions of control point grid.
  int m_Dims[3];

  /// Domain of control point grid in world coordinates.
  Types::Coordinate Domain[3];

  /// Array of spacings between the control points.
  Types::Coordinate Spacing[3];

  /// Array of spacings between the control points.
  Vector3D m_Offset;

  /// Get global scaling factor.
  virtual Types::Coordinate GetGlobalScaling() const
  { return this->GlobalScaling; }

  /// Initial affine transformation.
  cmtkGetSetMacro(AffineXform::SmartPtr,InitialAffineXform);

  /** Replace initial affine transformation.
   * If the new transformation is not given (or a NULL pointer), then
   * the new initial affine transformation is the identity transformation.
   * This means that the transformation is converted into its pure
   * nonrigid form.
   */
  void ReplaceInitialAffine( const AffineXform* newAffineXform = NULL );

  /// Concat affine transformation.
  void ConcatAffine( const AffineXform* affineXform );

  /// Get number of control points.
  size_t GetNumberOfControlPoints() const
  {
    return this->NumberOfControlPoints;
  }

protected:
  /// Number of control points.
  size_t NumberOfControlPoints;

  /** Inverted spacings between the control points.
   * These values are used for multiplication instead of division by those in
   * Spacing[].
   */
  Types::Coordinate InverseSpacing[3];

  /// Number of edge planes in the control point grid to keep unmoved.
  cmtkGetSetMacro(unsigned int,IgnoreEdge);

  /// Flag for fast but inaccurate computation.
  cmtkGetSetMacroDefault(bool,FastMode,true);

  /// Precomputed global scaling of initial affine transformation.
  Types::Coordinate GlobalScaling;

  /// Stored scale factors of the initial affine transformation.
  Types::Coordinate InverseAffineScaling[3];

  /// Offset of next control grid column.
  int nextI;

  /// Offset of next control grid row.
  int nextJ;

  /// Offset for next row and column.
  int nextIJ;

  /// Offset for next plane.
  int nextK;

  /// Offset for next plane and column.
  int nextIK;

  /// Offset for next plane and row.
  int nextJK;

  /// Offset for next plane, row, and column.
  int nextIJK;

public:
  /// Default constructor.
  WarpXform () : 
    m_InitialAffineXform( NULL ), 
    m_ActiveFlags( NULL )
  { 
    this->m_IgnoreEdge = 0; 
    this->m_FastMode = false; 
    IncompressibilityMap = DataGrid::SmartPtr( NULL );
    this->m_Dims[0] = this->m_Dims[1] = this->m_Dims[2] = 0;
    this->InverseSpacing[0] = this->InverseSpacing[1] = this->InverseSpacing[2] = 0.0;
  }

  /// Destructor.
  virtual ~WarpXform () {}

  /// Initialized internal data structures for new control point grid.
  virtual void InitGrid( const Types::Coordinate domain[3], const int dims[3] );

  /// Check whether coordinate is in domain of transformation.
  virtual bool InDomain( const Vector3D& v ) const 
  {
    return 
      ( v.XYZ[0] >= 0 ) && ( v.XYZ[0] <= Domain[0] ) &&
      ( v.XYZ[1] >= 0 ) && ( v.XYZ[1] <= Domain[1] ) &&
      ( v.XYZ[2] >= 0 ) && ( v.XYZ[2] <= Domain[2] );
  }
  
  /// Update internal representation.
  virtual void Update( const bool exactDelta = false );

  /// Refine control point grid, but maintain transformation exactly.
  virtual void Refine() {}

  /// Return warped vector.
  virtual Vector3D Apply ( const Vector3D& v ) const 
  {
    Vector3D w(v);
    this->ApplyInPlace( w );
    return w;
  }
  
  /** Return derivative of registration error with respect to one parameter.
   */
  virtual void GetDerivativeLandmarksMSD( double& lowerMSD, double& upperMSD, const MatchedLandmarkList* ll, const unsigned int idx, const Types::Coordinate step );

  /** Return inverse consistency.
   */
  virtual Types::Coordinate GetInverseConsistencyError( const WarpXform* inverse, const UniformVolume* volume, const Rect3D* voi = NULL ) const;
  
  /** Return derivative of inverse consistency.
   */
  virtual void GetDerivativeInverseConsistencyError
  ( double& lower, double& upper, const WarpXform* inverse, const UniformVolume* volume, const Rect3D* voi, 
    const unsigned int idx, const Types::Coordinate step );
  
  /// Set voxel-by-voxel map for incompressibility constraint.
  virtual void SetIncompressibilityMap( DataGrid::SmartPtr& incompressibility );

  /// Get the original position of a control point.
  virtual void GetOriginalControlPointPosition( Vector3D& v, const Types::Coordinate x, const Types::Coordinate y, const Types::Coordinate z) const 
  { 
    v.Set( this->m_Offset[0] + x*this->Spacing[0], this->m_Offset[1] + y*this->Spacing[1], this->m_Offset[2] + z*this->Spacing[2] );
  }
  
  /// Get the original position of a control point by index.
  virtual void GetOriginalControlPointPositionByOffset( Vector3D& v, const size_t offset ) const 
  {
    this->GetOriginalControlPointPosition( v, offset % this->m_Dims[0], (offset % (this->m_Dims[0]*this->m_Dims[1])) / this->m_Dims[0], offset / (this->m_Dims[0]*this->m_Dims[1]) ); 
  }

  /// Get shifted control point position.
  virtual void GetShiftedControlPointPosition( Vector3D& v, const int x, const int y, const int z ) const 
  { 
    this->GetShiftedControlPointPositionByOffset( v, x + this->m_Dims[0] * (y + this->m_Dims[1] * z ) );
  }

  /// Get shifted control point position by offset.
  virtual void GetShiftedControlPointPositionByOffset( Vector3D& v, const size_t offset ) const 
  { 
    v.Set( this->m_Parameters[offset*3], this->m_Parameters[offset*3+1], this->m_Parameters[offset*3+2] );
  }

  /// Set shifted control point position.
  virtual void SetShiftedControlPointPositionByOffset( const Vector3D& v, const int x, const int y, const int z ) const 
  { 
    this->SetShiftedControlPointPositionByOffset( v, x + this->m_Dims[0] * (y + this->m_Dims[1] * z ) );
  }

  /// Set shifted control point position by offset.
  virtual void SetShiftedControlPointPositionByOffset( const Vector3D& v, const size_t offset ) const 
  { 
    for ( int idx = 0; idx < 3; ++idx )
      this->m_Parameters[idx+offset*3] = v.XYZ[idx];
  }

  virtual Types::Coordinate GetParamStep( const size_t, const Types::Coordinate* volSize, const Types::Coordinate mmStep = 1 ) const;

  /// Free bitset for active parameter flags if it exists.
  void DeleteParameterActiveFlags();

  /// Set all parameters as active.
  void SetParameterActive();

  /** Set a particular parameter as active (or passive).
   *@param index Index of the parameter.
   *@param active Flag whether to set the parameter as active (non-zero) or
   * passive (zero).
   */
  void SetParameterActive( const size_t index, const bool active = true );

  /// Set a particular parameter as inactive.
  void SetParameterInactive( const size_t index );

  /// Set only those parameters as active that influence a given ROI.
  void SetParametersActive( const Rect3D& roi );
  
  /// Set parameters for one spatial dimension as active or inactive.
  void SetParametersActive( const int axis, const bool active = true );

  /** Set parameters for spatial dimensions as active.
   *@param axes This parameter defiend the activated dimensions in this
   * transformation. If it contains the characters x, y, or z, then the x, y,
   * and z-directions, respectively, are activated. The remaining directions
   * (if any) are deativated. The order of the axes is not relevant. Both 
   * upper- and lowercase characters will be accepted.
   */
  void SetParametersActive( const char* axes );

  /** Test whether a particular parameter is active.
   *@param index Index of the parameter.
   *@return Non-zero if and only if the given parameter is active.
   */
  int GetParameterActive( const size_t index ) const;

private:
  /** Flags for active (and passive) parameters.
   * This bitset contains one bit for each parameter in this transformation.
   * Every parameter is considered an active (variable) of passive (fixed)
   * parameter. Passive parameters are not considered for gradient computations
   * etc. and can therefore save a significant amount of computation time.
   */
  cmtkGetSetMacro(BitVector::SmartPtr,ActiveFlags);

  /** Voxel-by-voxel map of incompressibility constraint weight.
   */
  DataGrid::SmartPtr IncompressibilityMap;

  /// Friend declaration.
  friend class SplineWarpXformUniformVolume;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkWarpXform_h_included_
