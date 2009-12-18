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

#ifndef __cmtkImagePairNonrigidRegistrationFunctional_h_included_
#define __cmtkImagePairNonrigidRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkImagePairNonrigidRegistrationFunctional.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Parallel elastic registration functional.
 * This class provides multi-threaded implementations for the most 
 * time-consuming tasks performed by ImagePairNonrigidRegistrationFunctional and its
 * derived classes.
 */
template<class VM, class W> 
class ImagePairNonrigidRegistrationFunctional_Template
  /// Inherit from general image pair registration functional.
  : public ImagePairNonrigidRegistrationFunctional 
{
protected:
  /// Array of warp transformation objects for the parallel threads.
  SmartPointer<W> *m_ThreadWarp;

  /// Array of storage for simultaneously retrieving multiple deformed vectors.
  Vector3D **m_ThreadVectorCache;

  /** Number of actual parallel threads used for computations.
   * All duplicated data structures are generated with the multiplicity given
   * by this value. It is determined from Threads when the object is first
   * instanced. It cannot be changed afterwards.
   */
  size_t m_NumberOfThreads;

  /// Number of parallel tasks.
  size_t m_NumberOfTasks;

  /** Ground transformed volume.
   */
  typename Types::DataItem *m_WarpedVolume;

  /// Flag for forcing pixel values outside the floating image.
  bool m_ForceOutsideFlag;

  /// Rescaled byte value for forcing pixel values outside the floating image.
  typename Types::DataItem m_ForceOutsideValueRescaled;

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

public:
  /// This class.
  typedef ImagePairNonrigidRegistrationFunctional_Template<VM,W> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef ImagePairNonrigidRegistrationFunctional Superclass;

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

  /// Set flag and value for forcing values outside the floating image.
  virtual void SetForceOutside
  ( const bool flag = true, const Types::DataItem value = 0 )
  {
    this->m_ForceOutsideFlag = flag;
    this->m_ForceOutsideValueRescaled = this->Metric->DataY.ValueToIndex( value );
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

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v ) 
  {
    Warp->GetParamVector( v );
  }

  /// Constructor.
  ImagePairNonrigidRegistrationFunctional_Template<VM,W>( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating, const Interpolators::InterpolationEnum interpolation )
  {
    this->m_NumberOfThreads = ThreadPool::GlobalThreadPool.GetNumberOfThreads();
    this->m_NumberOfTasks = 4 * this->m_NumberOfThreads - 3;
    
    this->m_ThreadWarp = Memory::AllocateArray<typename W::SmartPtr>( this->m_NumberOfThreads );
    
    this->m_InfoTaskGradient.resize( this->m_NumberOfTasks );
    this->m_InfoTaskComplete.resize( this->m_NumberOfTasks );
    
    this->m_TaskMetric = Memory::AllocateArray<VM*>( this->m_NumberOfThreads );
    for ( size_t task = 0; task < this->m_NumberOfThreads; ++task )
      this->m_TaskMetric[task] = new VM( *(this->Metric) );
    
    this->m_ThreadVectorCache = Memory::AllocateArray<Vector3D*>( this->m_NumberOfThreads );
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      this->m_ThreadVectorCache[thread] = Memory::AllocateArray<Vector3D>( this->m_ReferenceDims[0] );
  }

  /** Destructor.
   * Free all per-thread data structures.
   */
  virtual ~ImagePairNonrigidRegistrationFunctional_Template<VM,W>() 
  {
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      if ( this->m_ThreadVectorCache[thread] ) 
	Memory::DeleteArray( this->m_ThreadVectorCache[thread] );
    Memory::DeleteArray( this->m_ThreadVectorCache );
    
    for ( size_t task = 0; task < this->m_NumberOfThreads; ++task )
      delete this->m_TaskMetric[task];
    Memory::DeleteArray( this->m_TaskMetric );
    
    Memory::DeleteArray( this->m_ThreadWarp );
  }

  /** Set warp transformation.
   * In the multi-threaded implementation, Warp[0] will be linked directly to
   * the given warp, while for all other threads a copy of the original object
   * is created by a call to WarpXform::Clone().
   */
  virtual void SetWarpXform ( WarpXform::SmartPtr& warp );
  
  /** Evaluate functional for the complete image data.
   * This function builds the pre-computed deformed floating image that is 
   * later used for rapid gradient computation.
   */
  typename Self::ReturnType EvaluateComplete ( CoordinateVector& v ) 
  {
    this->Metric->Reset();
    if ( ! this->WarpedVolume ) 
      this->WarpedVolume = Memory::AllocateArray<Types::DataItem>(  this->DimsX * this->DimsY * this->DimsZ  );
    
    this->m_ThreadWarp[0]->SetParamVector( v );

    const size_t numberOfTasks = std::min<size_t>( this->m_NumberOfTasks, this->DimsY * this->DimsZ );
    for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
      {
      InfoTaskComplete[taskIdx].thisObject = this;
      }

    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->TaskMetric[taskIdx]->Reset();
      }
    
    ThreadPool::GlobalThreadPool.Run( EvaluateCompleteThread, this->InfoTaskComplete );
    
    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->Metric->AddMetric( *(this->TaskMetric[taskIdx]) );
      }
    
    return this->WeightedTotal( this->Metric->Get(), this->m_ThreadWarp[0] );
  }

  /** Evaluate functional after change of a single parameter.
   *@param warp The current deformation.
   *@param metric The metric computed for the base-deformed volume.
   *@param voi Volume-of-Influence for the parameter under consideration.
   *@return The metric after recomputation over the given volume-of-influence.
   */
  typename Self::ReturnType EvaluateIncremental( const W* warp, VM *const localMetric, const Rect3D* voi, Vector3D *const vectorCache ) 
  {
    Vector3D *pVec;
    int pX, pY, pZ, offset, r;
    int fltIdx[3];
    Types::Coordinate fltFrac[3];

    int endLineIncrement = ( voi->startX + ( this->DimsX - voi->endX) );
    int endPlaneIncrement = this->DimsX * ( voi->startY + (this->DimsY - voi->endY) );
    
    const Types::DataItem unsetY = this->Metric->DataY.padding();
    localMetric->CopyUnsafe( *this->Metric );
    r = voi->startX + this->DimsX * ( voi->startY + this->DimsY * voi->startZ );
    for ( pZ = voi->startZ; pZ<voi->endZ; ++pZ ) 
      {
      for ( pY = voi->startY; pY<voi->endY; ++pY ) 
	{
	pVec = vectorCache;
	warp->GetTransformedGridSequence( pVec,voi->endX-voi->startX, voi->startX, pY, pZ );
	for ( pX = voi->startX; pX<voi->endX; ++pX, ++r, ++pVec ) 
	  {
	  // Remove this sample from incremental metric according to "ground warp" image.
	  const Types::DataItem sampleX = this->Metric->GetSampleX( r );
	  if ( this->WarpedVolume[r] != unsetY )
	    localMetric->Decrement( sampleX, this->WarpedVolume[r] );
	  
	  // Tell us whether the current location is still within the floating volume and get the respective voxel.
	  Vector3D::CoordMultInPlace( *pVec, this->FloatingInverseDelta );
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
  virtual typename Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const typename Self::ParameterType step = 1 ) 
  {
    const typename Self::ReturnType current = this->EvaluateComplete( v );

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

    ThreadPool::GlobalThreadPool.Run( EvaluateGradientThread, InfoTaskGradient );
    
    return current;
  }

  /// Evaluate functional.
  virtual typename Self::ReturnType EvaluateAt ( CoordinateVector& v )
  {
    this->m_ThreadWarp[0]->SetParamVector( v );
    return this->Evaluate();
  }

  virtual typename Self::ReturnType Evaluate ()
  {
    this->Metric->Reset();

    const size_t numberOfTasks = std::min<size_t>( this->m_NumberOfTasks, this->DimsY * this->DimsZ );
    for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
      {
      InfoTaskComplete[taskIdx].thisObject = this;
      }
    
    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->TaskMetric[taskIdx]->Reset();
      }
    
    ThreadPool::GlobalThreadPool.Run( EvaluateCompleteThread, this->InfoTaskComplete );
    
    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->Metric->AddMetric( *(this->TaskMetric[taskIdx]) );
      }
    
    return this->WeightedTotal( this->Metric->Get(), this->m_ThreadWarp[0] );
  }

private:
  /** Metric object for threadwise computation.
   * The objects in this array are the per-thread equivalent of the
   * ImagePairNonrigidRegistrationFunctional::IncrementalMetric object.
   */
  VM** TaskMetric;
  
  /// Consistency histogram objects for threadwise computation.
  JointHistogram<unsigned int>** ThreadConsistencyHistogram;
  
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

    SmartPointer<W>& myWarp = me->m_ThreadWarp[threadIdx];
    myWarp->SetParamVector( *info->Parameters );
    
    VM* threadMetric = me->TaskMetric[threadIdx];
    Vector3D *vectorCache = me->ThreadVectorCache[threadIdx];
    Types::Coordinate *p = myWarp->m_Parameters;
    
    Types::Coordinate pOld;
    double upper, lower;

    const Rect3D *voi = me->VolumeOfInfluence + taskIdx;
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
	upper = me->EvaluateIncremental( myWarp, threadMetric, voi, vectorCache );
	p[dim] = pOld - thisStep;
	lower = me->EvaluateIncremental( myWarp, threadMetric, voi, vectorCache );
	
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
    const W *warp = me->m_ThreadWarp[0];
    VM* threadMetric = me->TaskMetric[threadIdx];
    Vector3D *vectorCache = me->ThreadVectorCache[threadIdx];
    
    typename Types::DataItem* warpedVolume = me->WarpedVolume;
    const Types::DataItem unsetY = me->Metric->DataY.padding();
    
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
    
    int offset, r = rowFrom * me->DimsX;
    for ( pZ = pZfrom; (pZ < me->DimsZ) && rowsToDo; ++pZ ) 
      {
      for ( pY = pYfrom; (pY < me->DimsY) && rowsToDo; pYfrom = 0, ++pY, --rowsToDo ) 
	{
	warp->GetTransformedGridSequence( vectorCache, me->DimsX, 0, pY, pZ );
	pVec = vectorCache;
	for ( pX = 0; pX<me->DimsX; ++pX, ++r, ++pVec ) 
	  {
	  // Tell us whether the current location is still within the 
	  // floating volume and get the respective voxel.
	  Vector3D::CoordMultInPlace( *pVec, me->FloatingInverseDelta );
	  if ( me->FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	    {
	    // Compute data index of the floating voxel in the floating 
	    // volume.
	    offset = fltIdx[0] + me->FltDimsX * ( fltIdx[1] + me->FltDimsY*fltIdx[2] );
	    
	    // Continue metric computation.
	    warpedVolume[r] = me->Metric->GetSampleY(offset, fltFrac );
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
};

//@}

} // namespace cmtk

#endif // __cmtkImagePairNonrigidRegistrationFunctional_h_included_
