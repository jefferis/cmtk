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

#ifndef __cmtkVoxelMatchingElasticFunctional_h_included_
#define __cmtkVoxelMatchingElasticFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkVoxelMatchingFunctional.h>

#include <Base/cmtkVector.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkJointHistogram.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkMacros.h>
#include <Base/cmtkMathUtil.h>

#include <System/cmtkException.h>
#include <System/cmtkThreads.h>
#include <System/cmtkThreadPool.h>

#include <cassert>
#include <vector>

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
  virtual void SetWarpXform( SplineWarpXform::SmartPtr& warp ) = 0;

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

  /** Active coordinate directions.
   */
  cmtkGetSetMacroString(ActiveCoordinates);

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
  JointHistogram<unsigned int>::SmartPtr m_ConsistencyHistogram;

  /// Dimension of warp parameter vector
  size_t Dim;

  /** Parameter scaling vector.
   * This array holds the scaling factors for all warp parameters as returned
   * by the transformation class. These factors can be used to equalized all
   * parameter modifications during gradient computation etc.
   */
  std::vector<Types::Coordinate> StepScaleVector;

  /** Volume of influence table.
   * This array holds the precomputed volumes of influence for all 
   * transformation parameters. Six successive numbers per parameter define the
   * voxel range with respect to the reference colume grid that is affected by
   * the respective parameter.
   */
  DataGrid::RegionType *VolumeOfInfluence;

  /// Reference volume coordinate domain.
  UniformVolume::CoordinateRegionType m_ReferenceDomain;

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
  void SetInverseTransformation( typename W::SmartPtr& inverseTransformation ) 
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
  void WeightedDerivative( double& lower, double& upper, W& warp, const int param, const Types::Coordinate step ) const;

public:
  /// Get parameter stepping in milimeters.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    return Warp->GetParamStep( idx, Vector3D( this->FloatingSize ), mmStep );
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
  virtual void SetWarpXform ( typename W::SmartPtr& warp );

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
template<class VM> 
class VoxelMatchingElasticFunctional_Template 
  : public VoxelMatchingFunctional_Template<VM>,
    public VoxelMatchingElasticFunctional_WarpTemplate<SplineWarpXform> 
{
public:
  /// This class.
  typedef VoxelMatchingElasticFunctional_Template<VM> Self;

  /// Superclass.
  typedef VoxelMatchingElasticFunctional_WarpTemplate<SplineWarpXform> Superclass;

  /** Constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   *\param reference The reference (i.e. static) volume.
   *\param floating The floating (i.e. transformed) volume.
   */
  VoxelMatchingElasticFunctional_Template( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
    : VoxelMatchingFunctional_Template<VM>( reference, floating ), 
      VoxelMatchingElasticFunctional_WarpTemplate<SplineWarpXform>( reference, floating ),
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

    ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
    this->m_NumberOfThreads = threadPool.GetNumberOfThreads();
    this->m_NumberOfTasks = 4 * this->m_NumberOfThreads - 3;
    
    ThreadWarp.resize( this->m_NumberOfThreads );
    
    this->InfoTaskGradient.resize( this->m_NumberOfTasks );
    this->InfoTaskComplete.resize( this->m_NumberOfTasks );
    
    this->TaskMetric.resize( this->m_NumberOfThreads );
    for ( size_t task = 0; task < this->m_NumberOfThreads; ++task )
      this->TaskMetric[task] = new VM( *(this->Metric) );
    
    this->ThreadVectorCache = Memory::ArrayC::Allocate<Vector3D*>( this->m_NumberOfThreads );
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      this->ThreadVectorCache[thread] = Memory::ArrayC::Allocate<Vector3D>( this->ReferenceDims[0] );
  }

  /// Virtual destructor.
  virtual ~VoxelMatchingElasticFunctional_Template() 
  {
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      if ( ThreadVectorCache[thread] ) 
	Memory::ArrayC::Delete( this->ThreadVectorCache[thread] );
    Memory::ArrayC::Delete( this->ThreadVectorCache );
    
    for ( size_t task = 0; task < this->m_NumberOfThreads; ++task )
      delete this->TaskMetric[task];

    if ( WarpedVolume ) 
      Memory::ArrayC::Delete( WarpedVolume );
  }

  /// Set flag and value for forcing values outside the floating image.
  virtual void SetForceOutside
  ( const bool flag = true, const Types::DataItem value = 0 )
  {
    this->m_ForceOutsideFlag = flag;
    this->m_ForceOutsideValueRescaled = this->Metric->DataY.ValueToIndex( value );
  }

  /** Set warp transformation.
   * In the multi-threaded implementation, Warp[0] will be linked directly to
   * the given warp, while for all other threads a copy of the original object
   * is created by a call to WarpXform::Clone().
   */
  virtual void SetWarpXform ( SplineWarpXform::SmartPtr& warp ) 
  {
    this->Superclass::SetWarpXform( warp );
    
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread ) 
      {
      if ( this->Warp ) 
	{
	if ( thread ) 
	  {
	  ThreadWarp[thread] = SplineWarpXform::SmartPtr( this->Warp->Clone() );
	  ThreadWarp[thread]->RegisterVolume( this->ReferenceGrid );
	  } 
	else 
	  {
	  ThreadWarp[thread] = this->Warp;
	  }
	} 
      else
	{
	ThreadWarp[thread] = SplineWarpXform::SmartPtr::Null();
	}
      }
  }
  
  /** Evaluate functional after change of a single parameter.
   *\param warp The current deformation.
   *\param localMetric The local working metric.
   *\param voi Volume-of-Influence for the parameter under consideration.
   *\param vectorCache Pre-allocated storage for holding transformed vectors.
   *\return The metric after recomputation over the given volume-of-influence.
   */
  typename Self::ReturnType EvaluateIncremental( const SplineWarpXform& warp, VM *const localMetric, const DataGrid::RegionType& voi, Vector3D *const vectorCache ) 
  {
    Vector3D *pVec;
    int pX, pY, pZ, offset, r;
    int fltIdx[3];
    Types::Coordinate fltFrac[3];

    int endLineIncrement = ( voi.From()[0] + ( this->DimsX - voi.To()[0]) );
    int endPlaneIncrement = this->DimsX * ( voi.From()[1] + (this->DimsY - voi.To()[1]) );
    
    const typename VM::Exchange unsetY = this->Metric->DataY.padding();
    *localMetric = *this->Metric;
    r = voi.From()[0] + this->DimsX * ( voi.From()[1] + this->DimsY * voi.From()[2] );
    for ( pZ = voi.From()[2]; pZ<voi.To()[2]; ++pZ ) 
      {
      for ( pY = voi.From()[1]; pY<voi.To()[1]; ++pY ) 
	{
	pVec = vectorCache;
	warp.GetTransformedGridRow( voi.To()[0]-voi.From()[0], pVec, voi.From()[0], pY, pZ );
	for ( pX = voi.From()[0]; pX<voi.To()[0]; ++pX, ++r, ++pVec ) 
	  {
	  // Remove this sample from incremental metric according to "ground warp" image.
	  const typename VM::Exchange sampleX = this->Metric->GetSampleX( r );
	  if ( this->WarpedVolume[r] != unsetY )
	    localMetric->Decrement( sampleX, this->WarpedVolume[r] );
	  
	  // Tell us whether the current location is still within the floating volume and get the respective voxel.
	  *pVec *= this->FloatingInverseDelta;
	  if ( this->FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	    {
	    // Compute data index of the floating voxel in the floating volume.
	    offset = fltIdx[0] + this->FltDimsX * ( fltIdx[1] + this->FltDimsY * fltIdx[2] );
	    
	    // Continue metric computation.
	    localMetric->Increment( sampleX, this->Metric->GetSampleY(offset, fltFrac ) );
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
  
  /** Using OpenMP, update set of active and passive parameters.
   * This function computes local entropies in the neighborhood of all control
   * points of the Warp transformation. Those control points for which both
   * reference and floating image have less than half the maximum entropy in
   * this neighborhood as compared to the rest of the image are set passive.
   * The passive parameters are not considered for gradient computation and
   * therefore save significant computation time.
   */
  virtual void UpdateWarpFixedParameters();

  /// Compute functional value and gradient.
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const typename Self::ParameterType step = 1 ) 
  {
    const typename Self::ReturnType current = this->EvaluateAt( v );

    if ( this->m_AdaptiveFixParameters && this->WarpNeedsFixUpdate ) 
      {
      this->UpdateWarpFixedParameters();
      }
    
    // Make sure we don't create more threads than we have parameters.
    // Actually, we shouldn't create more than the number of ACTIVE parameters.
    // May add this at some point. Anyway, unless we have A LOT of processors,
    // we shouldn't really ever have more threads than active parameters :))
    const size_t numberOfTasks = std::min<size_t>( this->m_NumberOfTasks, this->Dim );
    
    for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
      {
      InfoTaskGradient[taskIdx].thisObject = this;
      InfoTaskGradient[taskIdx].Step = step;
      InfoTaskGradient[taskIdx].Gradient = g.Elements;
      InfoTaskGradient[taskIdx].BaseValue = current;
      InfoTaskGradient[taskIdx].Parameters = &v;
      }

    ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
    threadPool.Run( EvaluateGradientThread, InfoTaskGradient, numberOfTasks );
    
    return current;
  }

  /// Evaluate functional.
  virtual typename Self::ReturnType EvaluateAt ( CoordinateVector& v )
  {
    ThreadWarp[0]->SetParamVector( v );
    return this->Evaluate();
  }

  virtual typename Self::ReturnType Evaluate ()
  {
    this->Metric->Reset();
    if ( ! this->WarpedVolume ) 
      this->WarpedVolume = Memory::ArrayC::Allocate<typename VM::Exchange>(  this->DimsX * this->DimsY * this->DimsZ  );

    const size_t numberOfTasks = std::min<size_t>( this->m_NumberOfTasks, this->DimsY * this->DimsZ );
    for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
      {
      InfoTaskComplete[taskIdx].thisObject = this;
      }
    
    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->TaskMetric[taskIdx]->Reset();
      }
    
    ThreadPool::GetGlobalThreadPool().Run( EvaluateCompleteThread, this->InfoTaskComplete, numberOfTasks );
    
    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->Metric->AddMetric( *(this->TaskMetric[taskIdx]) );
      }
    
    return this->WeightedTotal( this->Metric->Get(), ThreadWarp[0] );
  }

private:
  /** Metric object for threadwise computation.
   * The objects in this array are the per-thread equivalent of the
   * VoxelMatchingElasticFunctional::IncrementalMetric object.
   */
  std::vector<VM*> TaskMetric;
  
#ifdef _OPENMP
  /// Consistency histogram objects for threadwise computation.
  std::vector<JointHistogram<unsigned int>::SmartPtr> m_ThreadConsistencyHistograms;
#endif // #ifdef _OPENMP
  
  /** Thread parameter block for incremental gradient computation.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  class EvaluateGradientTaskInfo 
  {
  public:
    /** Pointer to the functional object that created the thread. */
    Self *thisObject;
    /// Current parameter vector.
    CoordinateVector *Parameters;
    /// Current global coordinate stepping.
    typename Self::ParameterType Step;
    /// Pointer to gradient vector that is the target for computation results.
    Types::Coordinate *Gradient;
    /// Base functional value used for comparing new values to.
    double BaseValue;
  };
  
  /// Info blocks for parallel threads evaluating functional gradient.
  std::vector<typename Self::EvaluateGradientTaskInfo> InfoTaskGradient;
  
  /** Compute functional gradient as a thread.
   * This function (i.e., each thread) iterates over all parameters of the
   * current warp transformation. Among all active (i.e., not disabled)
   * parameters, it selects the ones that have an index with modulus
   * equal to the threads index when divided by the total number of threads.
   * For these parameters, the thread computes the partial derivative of the
   * functional by finite-difference approximation.
   */
  static void EvaluateGradientThread( void* arg, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t ) 
  {
    typename Self::EvaluateGradientTaskInfo *info = static_cast<typename Self::EvaluateGradientTaskInfo*>( arg );
    
    Self *me = info->thisObject;

    SplineWarpXform& myWarp = *(me->ThreadWarp[threadIdx]);
    myWarp.SetParamVector( *info->Parameters );
    
    VM* threadMetric = me->TaskMetric[threadIdx];
    Vector3D *vectorCache = me->ThreadVectorCache[threadIdx];
    Types::Coordinate *p = myWarp.m_Parameters;
    
    Types::Coordinate pOld;
    double upper, lower;

    const DataGrid::RegionType *voi = me->VolumeOfInfluence + taskIdx;
    for ( size_t dim = taskIdx; dim < me->Dim; dim+=taskCnt, voi+=taskCnt ) 
      {
      if ( me->StepScaleVector[dim] <= 0 ) 
	{
	info->Gradient[dim] = 0;
	}
      else
	{
	const typename Self::ParameterType thisStep = info->Step * me->StepScaleVector[dim];
	
	pOld = p[dim];
	
	p[dim] += thisStep;
	upper = me->EvaluateIncremental( myWarp, threadMetric, *voi, vectorCache );
	p[dim] = pOld - thisStep;
	lower = me->EvaluateIncremental( myWarp, threadMetric, *voi, vectorCache );
	
	p[dim] = pOld;
	me->WeightedDerivative( lower, upper, myWarp, dim, thisStep );
	
	if ( (upper > info->BaseValue ) || (lower > info->BaseValue) ) 
	  {
	  // strictly mathematically speaking, we should divide here by step*StepScaleVector[dim], but StepScaleVector[idx] is either zero or a constant independent of idx
	  info->Gradient[dim] = upper - lower;
	  } 
	else
	  {
	  info->Gradient[dim] = 0;
	  }
	}
      }
  }
  
  /** Thread parameter block for complete functional evaluation.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  class EvaluateCompleteTaskInfo 
  {
  public:
    /** Pointer to the functional object that created the thread. */
    Self *thisObject;
  };
  
  /** Info blocks for parallel threads evaluating complete functional. */
  std::vector<typename Self::EvaluateCompleteTaskInfo> InfoTaskComplete;
    
  /// Multi-threaded implementation of complete metric evaluation.
  static void EvaluateCompleteThread ( void *arg, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t ) 
  {
    typename Self::EvaluateCompleteTaskInfo *info = static_cast<typename Self::EvaluateCompleteTaskInfo*>( arg );
    
    Self *me = info->thisObject;
    const SplineWarpXform& warp = *(me->ThreadWarp[0]);
    VM* threadMetric = me->TaskMetric[threadIdx];
    Vector3D *vectorCache = me->ThreadVectorCache[threadIdx];
    
    typename VM::Exchange* warpedVolume = me->WarpedVolume;
    const typename VM::Exchange unsetY = me->Metric->DataY.padding();
    
    Vector3D *pVec;
    int pX, pY, pZ;
    
    int fltIdx[3];
    Types::Coordinate fltFrac[3];
    
    int rowCount = ( me->DimsY * me->DimsZ );
    int rowFrom = ( rowCount / taskCnt ) * taskIdx;
    int rowTo = ( taskIdx == (taskCnt-1) ) ? rowCount : ( rowCount / taskCnt ) * ( taskIdx + 1 );
    int rowsToDo = rowTo - rowFrom;
    
    int pYfrom = rowFrom % me->DimsY;
    int pZfrom = rowFrom / me->DimsY;
    
    int r = rowFrom * me->DimsX;
    for ( pZ = pZfrom; (pZ < me->DimsZ) && rowsToDo; ++pZ ) 
      {
      for ( pY = pYfrom; (pY < me->DimsY) && rowsToDo; pYfrom = 0, ++pY, --rowsToDo ) 
	{
	warp.GetTransformedGridRow( me->DimsX, vectorCache, 0, pY, pZ );
	pVec = vectorCache;
	for ( pX = 0; pX<me->DimsX; ++pX, ++r, ++pVec ) 
	  {
	  // Tell us whether the current location is still within the floating volume and get the respective voxel.
	  *pVec *= me->FloatingInverseDelta;
	  if ( me->FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	    {
	    // Compute data index of the floating voxel in the floating 
	    // volume.
	    const size_t offset = fltIdx[0] + me->FltDimsX * ( fltIdx[1] + me->FltDimsY*fltIdx[2] );
	    
	    // Continue metric computation.
	    warpedVolume[r] = me->Metric->GetSampleY(offset, fltFrac );
	    threadMetric->Increment( me->Metric->GetSampleX(r), warpedVolume[r] );
	    } 
	  else 
	    {
	    if ( me->m_ForceOutsideFlag )
	      {
	      warpedVolume[r] = me->m_ForceOutsideValueRescaled;
	      threadMetric->Increment( me->Metric->GetSampleX(r), warpedVolume[r] );
	      }
	    else
	      {
	      warpedVolume[r] = unsetY;
	      }
	    }
	  }
	}
      }
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
   *\see #EvaluateIncremental
   */
   SmartPointer<VM> IncrementalMetric;
   
   /// Shortcut variables for x, y, z dimension of the reference image.
   DataGrid::IndexType::ValueType DimsX, DimsY, DimsZ;
   
   /// Shorcut variables for x and y dimensions of the floating image.
   DataGrid::IndexType::ValueType FltDimsX, FltDimsY;

  /// Array of warp transformation objects for the parallel threads.
   std::vector<SplineWarpXform::SmartPtr> ThreadWarp;
  
  /// Array of storage for simultaneously retrieving multiple deformed vectors.
   Vector3D **ThreadVectorCache;
  
  /** Number of actual parallel threads used for computations.
   * All duplicated data structures are generated with the multiplicity given
   * by this value. It is determined from Threads when the object is first
   * instanced. It cannot be changed afterwards.
   */
  size_t m_NumberOfThreads;
  
  /// Number of parallel tasks.
  size_t m_NumberOfTasks;
};

/** Create functional from matching template.
 * This constructor function returns a pointer to a newly created elastic
 * matching functional. The functional is created from the template
 * corresponding to the given parameters.
 *\return A pointer to the newly created functional of NULL if creation failed.
 *\param metric Index of the voxel similarity measure to be used.
 *\param refVolume Reference volume data.
 *\param fltVolume Floating volume data.
 * template.
 */
VoxelMatchingElasticFunctional* CreateElasticFunctional( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume );

} // namespace

#endif // __cmtkVoxelMatchingElasticFunctional_h_included_
