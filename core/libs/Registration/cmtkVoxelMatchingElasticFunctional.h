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

#ifndef __cmtkVoxelMatchingElasticFunctional_h_included_
#define __cmtkVoxelMatchingElasticFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkVoxelMatchingFunctional.h>

#include <cmtkVector.h>
#include <cmtkWarpXform.h>
#include <cmtkJointHistogram.h>
#include <cmtkUniformVolume.h>

#include <cmtkMacros.h>
#include <cmtkMathUtil.h>
#include <cmtkException.h>

#include <assert.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#ifdef HAVE_VALUES_H
#  include <values.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Common base class for all elastic registration functionals.
 * This class holds all members that are not related to the effective metric
 * and therefore need not be present in the derived template class.
 */
class VoxelMatchingElasticFunctional : 
  /// Inherit basic voxel matching functions.
  public VoxelMatchingFunctional 
{
public:
  /// This class.
  typedef VoxelMatchingElasticFunctional Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef VoxelMatchingFunctional Superclass;

  /** Set Warp transformation.
   * This virtual function will be overridden by the derived classes that add
   * the actual warp transformation as a template parameters. It serves as a
   * common access point to update the warp transformation after construction
   * of the functional.
   */
  virtual void SetWarpXform( WarpXform::SmartPtr& warp ) = 0;

  /// Set flag and value for forcing pixels outside the floating image.
  virtual void SetForceOutside( const bool flag = true, const Types::DataItem value = 0 ) = 0;

  /// Destructor.
  virtual ~VoxelMatchingElasticFunctional ();

protected:
  /// Constructor.
  VoxelMatchingElasticFunctional( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating );

  /** Set active and passive warp parameters adaptively.
   * If this flag is set, the functional will adaptively determine active and
   * passive parameters of the warp transformation prior to gradient 
   * computation.
   */
  cmtkGetSetMacroDefault(bool,AdaptiveFixParameters,true);

  /** Set threshold factor for selecting passive warp parameters adaptively.
   * If the flag AdaptiveFixParameters is set, this value determines the
   * threshold by which active vs. passive parameters are selected. All
   * control points are set to passive for which the local region entropy is
   * below this factor times sum of min and max region entropy. The default
   * value is 0.5.
   */
  cmtkGetSetMacro(double,AdaptiveFixThreshFactor);

  /** Weight of the Jacobian constraint relative to voxel similarity measure.
   * If this is zero, only the voxel-based similarity will be computed.
   */
  cmtkGetSetMacroDefault(double,JacobianConstraintWeight,0);

  /** Weight of the rigidity constraint relative to voxel similarity measure.
   */
  cmtkGetSetMacroDefault(double,RigidityConstraintWeight,0);

  /** Map of rigidity weights constraint relative to voxel similarity measure.
   */
  cmtkGetSetMacro(DataGrid::SmartPtr,RigidityConstraintMap);

  /** Spatial map of relative (tissue-specific) incompressibility constraint.
   */
  cmtkGetSetMacro(DataGrid::SmartPtr,IncompressibilityMap);

  /** Weight of the grid energy relative to voxel similarity measure.
   * If this is zero, only the voxel-based similarity will be computed. If
   * equal to one, only the grid energy will be computed.
   */
  cmtkGetSetMacroDefault(double,GridEnergyWeight,0);

  /** Regularize the deformation.
   */
  cmtkGetSetMacroDefault(bool,Regularize,false);

  /** Warp's fixed parameters need to be updated.
   * This flag is set when the warp transformation is set or modified. It
   * signals that the active and passive parameters of the transformation
   * will have to be updated before the next gradient computation.
   */
  bool WarpNeedsFixUpdate;

  /// Histogram used for consistency computation.
  JointHistogram<unsigned int>::SmartPtr ConsistencyHistogram;

  /// Dimension of warp parameter vector
  size_t Dim;

  /** Parameter scaling vector.
   * This array holds the scaling factors for all warp parameters as returned
   * by the transformation class. These factors can be used to equalized all
   * parameter modifications during gradient computation etc.
   */
  Types::Coordinate *StepScaleVector;

  /** Volume of influence table.
   * This array holds the precomputed volumes of influence for all 
   * transformation parameters. Six successive numbers per parameter define the
   * voxel range with respect to the reference colume grid that is affected by
   * the respective parameter.
   */
  Rect3D *VolumeOfInfluence;

  /// Coordinate of the beginning of the reference colume crop area.
  Vector3D ReferenceFrom;

  /// Coordinate of the end of the reference colume crop area.
  Vector3D ReferenceTo;

  /// Storage for simultaneously retrieving multiple deformed vectors.
  Vector3D *VectorCache;
};

/** Template class for elastic registration functional.
 * This class incorporates all deformation-specific parts of the registration
 * functional for voxel-based non-rigid registration. The deformation type
 * itself is defined by the template-parameter W.
 *\note
 * The way multithreading is implemented by this class is as follows: A
 * fixed number of threads (currently two) is created. Each thread is assigned
 * a unique number from 0 through numberOfThreads-1. It then iterates over all
 * parameters of the current warp transformation object, counting the variable
 * (i.e., non-fixed) parameters. Out of these, it computes the partial 
 * derivative of the image similarity measure for those parameters for which
 * the modulus of the active parameter index divided by the number of threads
 * equals its own thread id.
 *\par
 * There are several advantages to this approach over others. First, if one
 * were to have one thread compute the upper and one the lower neighbor in
 * parameter space for each parameter, that would require 2 times the number
 * of parameters many threads, resulting in severy computation overhead. 
 * Instead, with the approach implemented here, we only need the given number
 * of threads which will then work in parallel as much as possible.
 */
template<class W> 
class VoxelMatchingElasticFunctional_WarpTemplate :
  /// Inherit from non-template base class.
  public VoxelMatchingElasticFunctional 
{
public:
  /// This class.
  typedef VoxelMatchingElasticFunctional_WarpTemplate<W> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef VoxelMatchingElasticFunctional Superclass;

  /// Pointer to the local warp transformation.
  typename W::SmartPtr Warp;

protected:
  /// Optional inverse transformation for inverse-consistent deformation.
  typename W::SmartPtr InverseTransformation;

  /// Weight for inverse consistency constraint.
  double InverseConsistencyWeight;

public:
  /// Set inverse transformation.
  void SetInverseTransformation( WarpXform::SmartPtr& inverseTransformation ) 
  {
    this->InverseTransformation = W::SmartPtr::DynamicCastFrom( inverseTransformation );
  }

  /// Set inverse consistency weight
  void SetInverseConsistencyWeight( const double inverseConsistencyWeight ) 
  {
    this->InverseConsistencyWeight = inverseConsistencyWeight;
  }

protected:

  /// Return weighted combination of voxel similarity and grid energy.
	typename Self::ReturnType WeightedTotal( const typename Self::ReturnType metric, const W* warp ) const 
  {
    double result = metric;
    if ( this->m_JacobianConstraintWeight > 0 ) 
      {
      result -= this->m_JacobianConstraintWeight * warp->GetJacobianConstraint();
      } 
    
    if ( this->m_RigidityConstraintWeight > 0 ) 
      {
      if ( this->m_RigidityConstraintMap )
	{
	result -= this->m_RigidityConstraintWeight * warp->GetRigidityConstraint( this->m_RigidityConstraintMap );
	}
      else
	{
	result -= this->m_RigidityConstraintWeight * warp->GetRigidityConstraint();
	}
      } 
    
    if ( this->m_GridEnergyWeight > 0 ) 
      {
      result -= this->m_GridEnergyWeight * warp->GetGridEnergy();
      }
    
    if ( !finite( result ) ) 
      return -FLT_MAX;
    
    if ( this->m_MatchedLandmarkList ) 
      {
      result -= this->m_LandmarkErrorWeight * warp->GetLandmarksMSD( this->m_MatchedLandmarkList );
      }

    if ( InverseTransformation ) 
      {
      result -= this->InverseConsistencyWeight * warp->GetInverseConsistencyError( this->InverseTransformation, this->ReferenceGrid );
      }
    
    return static_cast<typename Self::ReturnType>( result );
  }
  
  /// Return weighted combination of similarity and grid energy derivatives.
  void WeightedDerivative( double& lower, double& upper, typename W::SmartPtr& warp, const int param, const Types::Coordinate step ) const;

public:
  /// Get parameter stepping in milimeters.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    return Warp->GetParamStep( idx, FloatingSize, mmStep );
  }

  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const 
  {
    return Warp->ParamVectorDim();
  }

  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const 
  {
    return Warp->VariableParamVectorDim();
  }

  /** Set warp transformation.
   * In the multi-threaded implementation, Warp[0] will be linked directly to
   * the given warp, while for all other threads a copy of the original object
   * is created by a call to WarpXform::Clone().
   */
  virtual void SetWarpXform ( WarpXform::SmartPtr& warp );

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v ) 
  {
    Warp->GetParamVector( v );
  }

protected:
  /// Constructor.
  VoxelMatchingElasticFunctional_WarpTemplate( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
    : VoxelMatchingElasticFunctional( reference, floating ),
      Warp( NULL ), InverseTransformation( NULL )
  {}
  
  /// Dummy virtual destructor.
  virtual ~VoxelMatchingElasticFunctional_WarpTemplate() {}
};

/** Functional that evaluates a voxel-based similarity measure.
 * This class defines the type of functional that is optimized during
 * voxel-based registration. It holds references to reference and floating data
 * and computes similarity as well as its gradient w.r.t. a given
 * transformation.
 *
 * The metric to be optimized is given by a template parameter, therefore 
 * allowing inlined code to be generated for efficient evaluation.
 */
template<class VM, class W> class VoxelMatchingElasticFunctional_Template 
  : public VoxelMatchingFunctional_Template<VM>,
    public VoxelMatchingElasticFunctional_WarpTemplate<W> 
{
public:
  /// This class.
  typedef VoxelMatchingElasticFunctional_Template<VM,W> Self;

  /** Constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   *@param reference The reference (i.e. static) volume.
   *@param floating The floating (i.e. transformed) volume.
   */
  VoxelMatchingElasticFunctional_Template( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
    : VoxelMatchingFunctional_Template<VM>( reference, floating ), 
      VoxelMatchingElasticFunctional_WarpTemplate<W>( reference, floating ),
      m_ForceOutsideFlag( false ),
      m_ForceOutsideValueRescaled( 0 )
  {
    IncrementalMetric = typename VM::SmartPtr( new VM( *(this->Metric) ) );

    WarpedVolume = NULL;

    DimsX = this->ReferenceGrid->GetDims()[0];
    DimsY = this->ReferenceGrid->GetDims()[1];
    DimsZ = this->ReferenceGrid->GetDims()[2];

    FltDimsX = this->FloatingGrid->GetDims()[0];
    FltDimsY = this->FloatingGrid->GetDims()[1];
  }

  /// Virtual destructor.
  virtual ~VoxelMatchingElasticFunctional_Template() 
  {
    if ( WarpedVolume ) delete[] WarpedVolume;
  }

  /// Set flag and value for forcing values outside the floating image.
  virtual void SetForceOutside
  ( const bool flag = true, const Types::DataItem value = 0 )
  {
    this->m_ForceOutsideFlag = flag;
    this->m_ForceOutsideValueRescaled = this->Metric->DataY.ValueToIndex( value );
  }

  /** Evaluate functional after change of a single parameter.
   *@param warp The current deformation.
   *@param metric The metric computed for the base-deformed volume.
   *@param voi Volume-of-Influence for the parameter under consideration.
   *@return The metric after recomputation over the given volume-of-influence.
   */
  typename Self::ReturnType EvaluateIncremental( const W* warp, SmartPointer<VM>& localMetric, const Rect3D* voi ) 
  {
    Vector3D *pVec;
    int pX, pY, pZ, offset, r;
    int fltIdx[3];
    Types::Coordinate fltFrac[3];

    int endLineIncrement = ( voi->startX + (DimsX - voi->endX) );
    int endPlaneIncrement = DimsX * ( voi->startY + (DimsY - voi->endY) );
    
    const typename VM::Exchange unsetY = this->Metric->DataY.padding();
    localMetric->CopyUnsafe( *this->Metric );
    r = voi->startX + DimsX * ( voi->startY + DimsY * voi->startZ );
    for ( pZ = voi->startZ; pZ<voi->endZ; ++pZ ) 
      {
      for ( pY = voi->startY; pY<voi->endY; ++pY ) 
	{
	pVec = this->VectorCache;
	warp->GetTransformedGridSequence( pVec, voi->endX-voi->startX, voi->startX, pY, pZ );
	for ( pX = voi->startX; pX<voi->endX; ++pX, ++r, ++pVec ) 
	  {
	  // Remove this sample from incremental metric according to "ground warp" image.
	  const typename VM::Exchange sampleX = this->Metric->GetSampleX( r );
	  if ( this->WarpedVolume[r] != unsetY )
	    localMetric->Decrement( sampleX, WarpedVolume[r] );
	    
	  // Tell us whether the current location is still within the floating volume and get the respective voxel.
	  Vector3D::CoordMultInPlace( *pVec, this->FloatingInverseDelta );
	  if ( this->FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	    {
	    // Compute data index of the floating voxel in the floating volume.
	    offset = fltIdx[0] + FltDimsX * ( fltIdx[1] + FltDimsY * fltIdx[2] );
	    
	    // Continue metric computation.
	    localMetric->Increment( sampleX, this->Metric->GetSampleY( offset, fltFrac ) );
	    } 
	  else 
	    {
	    if ( this->m_ForceOutsideFlag )
	      {
	      localMetric->Increment( sampleX, this->m_ForceOutsideValueRescaled );
	      }
	    }
	  }
	r += endLineIncrement;
	}
      r += endPlaneIncrement;
      }
    
    return localMetric->Get();
  }

  /** Update set of active and passive parameters.
   * This function computes local entropies in the neighborhood of all control
   * points of the Warp transformation. Those control points for which both
   * reference and floating image have less than half the maximum entropy in
   * this neighborhood as compared to the rest of the image are set passive.
   * The passive parameters are not considered for gradient computation and
   * therefore save significant computation time.
   */
  void UpdateWarpFixedParameters();

  /// Compute functional value and gradient.
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 ) 
  {
    const typename Self::ReturnType current = this->EvaluateAt( v );

    if ( this->m_AdaptiveFixParameters && this->WarpNeedsFixUpdate ) 
      {
      this->UpdateWarpFixedParameters();
      }
    
    Types::Coordinate *p = this->Warp->m_Parameters;

    Types::Coordinate pOld;
    double upper, lower;
    Rect3D *voi = this->VolumeOfInfluence;
    for ( size_t dim = 0; dim < this->Dim; ++dim, ++voi ) 
      {
      if ( this->StepScaleVector[dim] <= 0 ) 
	{
	g[dim] = 0;
	} 
      else
	{
	pOld = p[dim];
	
	Types::Coordinate thisStep = step * this->StepScaleVector[dim];
	
	p[dim] += thisStep;
	upper = EvaluateIncremental( this->Warp, IncrementalMetric, voi );
	
	p[dim] = pOld - thisStep;
	lower = EvaluateIncremental( this->Warp, IncrementalMetric, voi );
	
	p[dim] = pOld;
	this->WeightedDerivative( lower, upper, this->Warp, dim, thisStep );
	
	if ( (upper > current ) || (lower > current) ) 
	  {
	  g[dim] = upper-lower;
	  } 
	else 
	  {
	  g[dim] = 0;
	  }
	}
      }
    
    return current;
  }
  
  /// Evaluate functional.
  virtual typename Self::ReturnType EvaluateAt( CoordinateVector& v )
  {
    this->Warp->SetParamVector( v );
    return this->Evaluate();
  }

  /// Evaluate functional.
  virtual typename Self::ReturnType Evaluate()
  {
    if ( ! WarpedVolume ) 
      WarpedVolume = Memory::AllocateArray<typename VM::Exchange>(  DimsX * DimsY * DimsZ  );

    this->Metric->Reset();
    const typename VM::Exchange unsetY = this->Metric->DataY.padding();

    Vector3D *pVec;
    int pX, pY, pZ;
    int fltIdx[3];
    Types::Coordinate fltFrac[3];

    int r = 0;
    for ( pZ = 0; pZ<DimsZ; ++pZ ) 
      {
      for ( pY = 0; pY<DimsY; ++pY )
	{
	this->Warp->GetTransformedGridSequence( this->VectorCache, DimsX, 0, pY, pZ );
	pVec = this->VectorCache;
	for ( pX = 0; pX<DimsX; ++pX, ++r, ++pVec )
	  {
	  // Tell us whether the current location is still within the floating volume and get the respective voxel.
	  Vector3D::CoordMultInPlace( *pVec, this->FloatingInverseDelta );
	  if ( this->FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	    {
	    // Compute data index of the floating voxel in the
	    // floating volume.
	    const size_t offset = fltIdx[0] + FltDimsX * ( fltIdx[1] + FltDimsY * fltIdx[2] );
	    
	    // Continue metric computation.
	    this->WarpedVolume[r] = this->Metric->GetSampleY( offset, fltFrac );
	    this->Metric->Increment( this->Metric->GetSampleX( r ), this->WarpedVolume[r] );
	    }
	  else 
	    {
	    if ( this->m_ForceOutsideFlag )
	      {
	      this->WarpedVolume[r] = this->m_ForceOutsideValueRescaled;
	      this->Metric->Increment( this->Metric->GetSampleX( r ), this->WarpedVolume[r] );
	      }
	    else
	      {
	      this->WarpedVolume[r] = unsetY;
	      }
	    }
	  }
	}
      }
    
    return this->WeightedTotal( this->Metric->Get(), this->Warp );
  }

protected:
  /** Ground transformed volume.
   */
  typename VM::Exchange *WarpedVolume;

  /// Flag for forcing pixel values outside the floating image.
  bool m_ForceOutsideFlag;

  /// Rescaled byte value for forcing pixel values outside the floating image.
  typename VM::Exchange m_ForceOutsideValueRescaled;

  /** Metric object for incremental computation.
   * Before computing the incremental metric after change of one parameter,
   * the global metric is copied to this object. It is then used for in-place
   * application of all necessary changes, leaving the original metric intact.
   *@see #EvaluateIncremental
   */
   SmartPointer<VM> IncrementalMetric;
   
   /// Shortcut variables for x, y, z dimension of the reference image.
   GridIndexType DimsX, DimsY, DimsZ;
   
   /// Shorcut variables for x and y dimensions of the floating image.
   GridIndexType FltDimsX, FltDimsY;
};

/** Create functional from matching template.
 * This constructor function returns a pointer to a newly created elastic
 * matching functional. The functional is created from the template
 * corresponding to the given parameters.
 *@return A pointer to the newly created functional of NULL if creation failed.
 *@param metric Index of the voxel similarity measure to be used.
 *@param warp Type of elastic transformation: Linear (0) or Spline (1).
 *@param refVolume Reference volume data.
 *@param modVolume Floating volume data.
 * template.
 */
VoxelMatchingElasticFunctional* CreateElasticFunctional( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& modVolume );

} // namespace

#endif // __cmtkVoxelMatchingElasticFunctional_h_included_
