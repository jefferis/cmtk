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

#ifndef __cmtkSplineWarpXform_h_included_
#define __cmtkSplineWarpXform_h_included_

#include <cmtkconfig.h>

#include <cmtkWarpXform.h>

#include <cmtkVector.h>
#include <cmtkVector3D.h>
#include <cmtkRect3D.h>
#include <cmtkAffineXform.h>
#include <cmtkCubicSpline.h>

#include <assert.h>
#include <vector>
#include <algorithm>
  
#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** B-spline-based local deformation.
 */
class SplineWarpXform :
  /// Inherit generic warp interface and basic functions.
  public WarpXform 
{
public:
  /// This class.
  typedef SplineWarpXform Self;

  /// Parent class.
  typedef WarpXform Superclass;

  /// Smart pointer to SplineWarpXform.
  typedef SmartPointer<Self> SmartPtr;

  /** Construct empty new warp.
   * This is mostly useful for the Clone() function that can subsequently
   * copy all data structures into a newly created empty object. 
   */
  SplineWarpXform();

  /** Construct new warp from volume and control grid density.
   */
  SplineWarpXform( const UniformVolume *volume,  const Types::Coordinate delta, const bool exactDelta = false );

  /** Construct new warp from volume size and control grid density.
   */
  SplineWarpXform( const Types::Coordinate domain[3], const Types::Coordinate delta, const AffineXform *initialXform = NULL, const bool exactDelta = false );

  /** Initialize warp from volume size and control grid density.
   */
  void Init( const Types::Coordinate domain[3], const Types::Coordinate delta, const AffineXform *initialXform = NULL, const bool exactDelta = false );
  
  /** Construct new warp from volume size, grid dimensions and parameters
   */
  SplineWarpXform( const Types::Coordinate domain[3], const int dims[3], CoordinateVector::SmartPtr& parameters, const AffineXform *initialXform = NULL );

  /// Initialize control point positions, potentially with affine displacement.
  void InitControlPoints( const AffineXform* affineXform = NULL );

  /// Get the original position of a control point.
  virtual void GetOriginalControlPointPosition( Vector3D& v, const Types::Coordinate x, const Types::Coordinate y, const Types::Coordinate z) const 
  { 
    v.Set( (x-1)*this->Spacing[0], (y-1)*this->Spacing[1], (z-1)*this->Spacing[2] );
  }
  
  /// Get the original position of a control point by index.
  virtual void GetOriginalControlPointPositionByOffset( Vector3D& v, const size_t offset ) const 
  {
    this->GetOriginalControlPointPosition( v, offset % this->m_Dims[0], (offset % (this->m_Dims[0]*this->m_Dims[1])) / this->m_Dims[0], offset / (this->m_Dims[0]*this->m_Dims[1]) ); 
  }

  /// Update internal representation.
  virtual void Update( const bool exactDelta = false );

  /// Clone transformation.
  virtual SplineWarpXform* Clone () const;

  /** Create inverse transformation.
   * This function returns NULL as there is no explicit inverse of a spline
   * warp transformation.
   */
  virtual SplineWarpXform* MakeInverse() const 
  {
    return NULL;
  }
  
  /// Refine transformation.
  virtual void Refine ( const int factor = 2 );

  /// Return grid bending energy.
  virtual Types::Coordinate GetGridEnergy() const;

  /** Return grid bending energy at one control point.
   *@param cp The control point where the bending energy is to be evaluated.
   */
  virtual Types::Coordinate GetGridEnergy( const Types::Coordinate *cp ) const;

  /** Return grid bending energy at arbitrary location.
   */
  virtual Types::Coordinate GetGridEnergy( const Vector3D& v ) const;

  /// Return derivative of grid energy with respect to one parameter.
  virtual void GetGridEnergyDerivative( double& lower, double& upper, const int param, const Types::Coordinate step ) const;

  /// Compute Jacobian determinant at a certain location.
  virtual Types::Coordinate GetJacobianDeterminant ( const Vector3D& v ) const;

  /// Compute Jacobian determinant at a certain reference image pixel.
  virtual Types::Coordinate GetJacobianDeterminant ( const int x, const int y, const int z ) const;

  /// Compute sequence of Jacobian determinants from given grid location.
  virtual void GetJacobianDeterminantSequence( double *const values, const int x, const int y, const int z, const size_t numberOfPoints = 1 ) const;

  /** Compute Jacobian determinant at a certain control point.
   * As we evaluate the spline polynomials and their derivatives at a known
   * control point, they reduce to constants which allows for a particularly
   * efficient computation.
   *@param cp Pointer to the spline control point at which the Jacobian
   * determinant is to be evaluated. This pointer is assumed to be pointing
   * to the 0th-component (ie, x) of a non-peripheral control point. Callers
   * have to make sure that this is true.
   */
  Types::Coordinate JacobianDeterminant ( const Types::Coordinate *cp ) const;

  /// Return Jacobian constraint of the current transformation grid.
  virtual Types::Coordinate GetJacobianConstraint() const;

  /** Return Jacobian folding constraint of the current transformation grid. */
  virtual Types::Coordinate GetJacobianFoldingConstraint() const;

  /// Return rigidity constraint of the current transformation grid.
  virtual Types::Coordinate GetRigidityConstraint() const;

  /// Return rigidity constraint of the current transformation grid with local weights.
  virtual Types::Coordinate GetRigidityConstraint( const DataGrid* weightMap ) const;

  /** Return sparse Jacobian constraint of the current transformation grid.
   * Unlike GetJacobianConstraint(), this function does not evaluate the
   * Jacobian constraint at each voxel of the registered reference image; it
   * merely performs the evaluation at the control points of the deformation,
   * thereby greatly increasing computational performance at the cost of less
   * accurate and stable computation.
   */
  virtual Types::Coordinate GetJacobianConstraintSparse() const;

  /** Return sparse rigidity constraint of the current transformation grid.
   */
  virtual Types::Coordinate GetRigidityConstraintSparse() const;

  /// Return derivative of Jacobian constraint with respect to one parameter.
  virtual void GetJacobianConstraintDerivative( double& lower, double& upper, const int param, const Rect3D&, const Types::Coordinate step ) const;

  /// Return derivative of Jacobian constraint with respect to one parameter.
  virtual void GetJacobianFoldingConstraintDerivative( double& lower, double& upper, const int param, const Rect3D&, const Types::Coordinate step ) const;

  /// Return derivative of Jacobian constraint with respect to one parameter.
  virtual void GetJacobianConstraintDerivative( double& lower, double& upper, const int param, const Types::Coordinate step ) const;
  
  /// Return derivative of rigidity constraint with respect to one parameter.
  virtual void GetRigidityConstraintDerivative( double& lower, double& upper, const int param, const Rect3D&, const Types::Coordinate step ) const;
  
  /// Return derivative of rigidity constraint with respect to one parameter.
  virtual void GetRigidityConstraintDerivative( double& lower, double& upper, const int param, const Rect3D&, const Types::Coordinate step,
						const DataGrid* weightMap ) const;

  /// Return derivative of rigidity constraint with respect to one parameter.
  virtual void GetRigidityConstraintDerivative( double& lower, double& upper, const int param, const Types::Coordinate step ) const;

  /** Return inverse consistency.
   */
  virtual Types::Coordinate GetInverseConsistencyError( const WarpXform* inverse, const UniformVolume* volume, const Rect3D* voi = NULL ) const;
  
  /** Return origin of warped vector.
   * Note that since this class of transformation is not closed under inversion
   * this function computes only a more or less accurate numerical 
   * approximation to the actual origin of a warped vector. Note also that this
   * computation is everything but computationally efficient.
   *@return True is the given inverse was succesfully comuted, false if the
   * given warped vector was outside the target domain of this transformation.
   */
  virtual bool ApplyInverse ( const Vector3D& v, Vector3D& u, const Types::Coordinate accuracy = 0.01  ) const;

  /** Return origin of warped vector.
   * Note that since this class of transformation is not closed under inversion
   * this function computes only a more or less accurate numerical 
   * approximation to the actual origin of a warped vector. Note also that this
   * computation is everything but computationally efficient.
   *@return True is the given inverse was succesfully comuted, false if the
   * given warped vector was outside the target domain of this transformation.
   */
  virtual bool ApplyInverseInPlace( Vector3D& v, const Types::Coordinate accuracy = 0.01  ) const;

  /** Return origin of warped vector.
   * Note that since this class of transformation is not closed under inversion
   * this function computes only a more or less accurate numerical 
   * approximation to the actual origin of a warped vector. Note also that this
   * computation is everything but computationally efficient.
   *@param initial Initial estimate for the original location. Search goes
   * from here. This is useful for looking up the original locations of
   * a large number of closely located vectors, for example all pixels in an
   * image.
   *@return True is the given inverse was succesfully comuted, false if the
   * given warped vector was outside the target domain of this transformation.
   */
  virtual bool ApplyInverseInPlaceWithInitial( Vector3D& v, const Vector3D& initial, const Types::Coordinate accuracy = 0.01 ) const;

  /// Replace existing vector with transformed location.
  virtual void ApplyInPlace( Vector3D& v ) const 
  {
    Types::Coordinate r[3], f[3];
    int grid[3];
    
    // Do some precomputations.
    { 
    for ( int dim = 0; dim<3; ++dim ) 
      {
      // This is the (real-valued) index of the control point grid cell the
      // given location is in.
      r[dim] = this->InverseSpacing[dim] * v.XYZ[dim];
      // This is the actual cell index.
      grid[dim] = std::min<int>( static_cast<int>( r[dim] ), this->m_Dims[dim]-4 );
      // And here's the relative position within the cell.
      f[dim] = r[dim] - grid[dim];
      }
    }

    // Create a pointer to the front-lower-left corner of the c.p.g. cell.
    const Types::Coordinate* coeff = this->m_Parameters + 3 * ( grid[0] + this->m_Dims[0] * (grid[1] + this->m_Dims[1] * grid[2]) );

    for ( int dim = 0; dim<3; ++dim ) 
      {
      Types::Coordinate mm = 0;
      const Types::Coordinate *coeff_mm = coeff;
      
      // Loop over 4 c.p.g. planes in z-direction.
      for ( int m = 0; m < 4; ++m ) 
	{
	Types::Coordinate ll = 0;
	const Types::Coordinate *coeff_ll = coeff_mm;
	
	// Loop over 4 c.p.g. planes in y-direction.
	for ( int l = 0; l < 4; ++l ) 
	  {
	  Types::Coordinate kk = 0;
	  const Types::Coordinate *coeff_kk = coeff_ll;
	  
	  // Loop over 4 c.p.g. planes in x-direction.
	  for ( int k = 0; k < 4; ++k, coeff_kk+=3 ) 
	    {
	    kk += CubicSpline::ApproxSpline( k, f[0] ) * (*coeff_kk);
	    }
	  ll += CubicSpline::ApproxSpline( l, f[1] ) * kk;
	  coeff_ll += nextJ;
	  }	
	mm += CubicSpline::ApproxSpline( m, f[2] ) * ll;
	coeff_mm += nextK;
	}
      v.XYZ[dim] = mm;
      ++coeff;
      }
  }

  /// Get volume influenced by one parameter.
  virtual void GetVolumeOfInfluence( const size_t idx, const Vector3D&, const Vector3D&, Vector3D&, Vector3D&, const int = -1 ) const;
  
  /// Unegister axes points, ie free all internal data structures.
  virtual void UnRegisterVolume();
  
  /// Get a grid point from the deformed grid.
  virtual void GetTransformedGrid( Vector3D& v, const int idxX, const int idxY, const int idxZ ) const;
  
  /// Get a sequence of grid points from the deformed grid. 
  virtual void GetTransformedGridSequence( Vector3D *const v, const int numPoints, const int idxX, const int idxY, const int idxZ ) const;
  
  /// Apply transformation to all coordinate triples stored in an array.
  virtual void ApplyToAll( CoordinateVector& v, BitVector& valid, const bool inverse = false, const Types::Coordinate epsilon = 0.0, const int* gridDims = NULL );
  
  /// Get parameter stepping.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate* volSize, const Types::Coordinate mmStep = 1 ) const 
  {
    return 4 * WarpXform::GetParamStep( idx, volSize, mmStep );
  }
  
  /** Get the deformed position of a transformation control point.
   *\note This function does not return the shifted control point position,
   * but rather it applies the current transformation to the given control
   * point. It does so more efficiently than applying the transformation to
   * the explicit 3D coordinate of the control point, because most spline
   * coefficients vanish at control points and need not be considered.
   */
  virtual Vector3D& GetDeformedControlPointPosition( Vector3D&, const int, const int, const int ) const;
  
  /** Return array of pre deformation vectors.
   * The newly alocated data array contains the control point positions
   * after deformation without affine components.
   *\param includeScale If this flag is set (default: off), then the scale
   * components of the affine transformation remain in the deformation.
   */
  Types::Coordinate* GetPureDeformation( const bool includeScale = false ) const;

  /// Get local Jacobian.
  virtual CoordinateMatrix3x3 GetJacobian( const Vector3D& v ) const;

  /// Get local Jacobian into existing matrix.
  virtual void GetJacobian( const Vector3D& v, CoordinateMatrix3x3& J ) const;

  /// Get local Jacobian at control point into existing matrix.
  virtual void GetJacobianAtControlPoint( const Types::Coordinate* cp, CoordinateMatrix3x3& J ) const;

  /// Get sequence of Jacobians for pixel row.
  virtual void GetJacobianSequence( CoordinateMatrix3x3 *const array, const int x, const int y, const int z, const size_t numberOfPoints ) const;
  
  /// Get local strain tensor.
  virtual CoordinateMatrix3x3 GetStrainTensor( const Vector3D& v ) const;

  /// Get local strain tensor into existing matrix.
  virtual void GetStrainTensor( const Vector3D& v, CoordinateMatrix3x3& J ) const;

private:
  /// Register axes points of the volume to be deformed.
  void RegisterVolumePoints ( const int[3], const Types::Coordinate *const* );

  /// Register axes points of the volume to be deformed.
  void RegisterVolumePoints( const int[3], const Types::Coordinate[3], const Types::Coordinate[3] );

  /// Register a single axis of the uniform volume to be deformed.
  void RegisterVolumeAxis ( const int, const Types::Coordinate delta, const Types::Coordinate origin, const int, const Types::Coordinate, std::vector<int>& g, 
			    std::vector<Types::Coordinate>& spline, std::vector<Types::Coordinate>& dspline );

  /// Return rigidity constraint based on given Jacobian matrix.
  Types::Coordinate GetRigidityConstraint( const CoordinateMatrix3x3& J ) const;

protected:
  /// Dimensions of the volume image linked to this transformation.
  int VolumeDims[3];

  /**@name Precomputed grid indices.
   * These arrays hold the precomputed grid indices of the deformed grid's
   * voxels with respect to the control point grid of this deformation.
   */
  //@{
  /// x-axis.
  std::vector<int> gX;
  /// y-axis.
  std::vector<int> gY;
  /// z-axis.
  std::vector<int> gZ;
  //@}

  /**@name Precomputed spline coefficients.
   * These arrays hold the precomputed spline coefficients for deforming the
   * voxel locations in the associated deformed grid.
   */
  //@{
  /// x-axis.
  std::vector<Types::Coordinate> splineX;
  /// y-axis.
  std::vector<Types::Coordinate> splineY;
  /// z-axis.
  std::vector<Types::Coordinate> splineZ;
  //@}

  /**@name Precomputed derivative spline coefficients.
   * These arrays hold the precomputed derivatives of the spline coefficients.
   * This allows for rapid evaluation of the Jacobian determinant.
   */
  //@{
  /// x-axis.
  std::vector<Types::Coordinate> dsplineX;
  /// y-axis.
  std::vector<Types::Coordinate> dsplineY;
  /// z-axis.
  std::vector<Types::Coordinate> dsplineZ;
  //@}

  /// Relative offsets of all control points in a 4 x 4 x 4 neighborhood.
  int GridPointOffset[48];

  /** Initialize internal data structures.
   * This function is called from the various destructors to avoid unnecessary
   * duplication of code.
   */
  void Init ();

  /** Thread parameter block for volume resampling.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  typedef struct 
  {
    /// Pointer to the functional object that created the thread.
    const SplineWarpXform *thisObject;
    /// Unique index of this thread instance among all threads.
    int ThisThreadIndex;
    /// Total number of threads created.
    int NumberOfThreads;
    /// Constraint for subvolume handled by this thread.
    Types::Coordinate Constraint;
  } JacobianConstraintThreadInfo;
  
  static CMTK_THREAD_RETURN_TYPE GetJacobianConstraintThreads( void *arg );

  static CMTK_THREAD_RETURN_TYPE GetJacobianFoldingConstraintThreads( void *arg );

  /// Find nearest (after deformation) control point.
  void FindClosestControlPoint( const Vector3D& v, Vector3D& cp ) const;

  /// Friend declaration.
  friend class SplineWarpXformUniformVolume;
};

} // namespace

#endif // #ifndef __cmtkSplineWarpXform_h_included_
