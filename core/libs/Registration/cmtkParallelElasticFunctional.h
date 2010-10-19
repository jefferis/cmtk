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

#ifndef __cmtkParallelElasticFunctional_h_included_
#define __cmtkParallelElasticFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkVoxelMatchingElasticFunctional.h>

#include <System/cmtkThreads.h>
#include <System/cmtkThreadPool.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Parallel elastic registration functional.
 * This class provides multi-threaded implementations for the most 
 * time-consuming tasks performed by VoxelMatchingElasticFunctional and its
 * derived classes.
 */
template<class VM> 
class ParallelElasticFunctional
  /// Inherit from non-parallel functional.
  : public VoxelMatchingElasticFunctional_Template<VM> 
{
protected:
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

public:
  /// This class.
  typedef ParallelElasticFunctional<VM> Self;

  /// Superclass.
  typedef VoxelMatchingElasticFunctional_Template<VM> Superclass;

  /// Constructor.
  ParallelElasticFunctional ( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating ) :
    VoxelMatchingElasticFunctional_Template<VM>( reference, floating )
  {
    ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
    this->m_NumberOfThreads = threadPool.GetNumberOfThreads();
    this->m_NumberOfTasks = 4 * this->m_NumberOfThreads - 3;
    
    ThreadWarp.resize( this->m_NumberOfThreads );
    
    this->InfoTaskGradient.resize( this->m_NumberOfTasks );
    this->InfoTaskComplete.resize( this->m_NumberOfTasks );
    
    this->TaskMetric = Memory::AllocateArray<VM*>( this->m_NumberOfThreads );
    for ( size_t task = 0; task < this->m_NumberOfThreads; ++task )
      this->TaskMetric[task] = new VM( *(this->Metric) );
    
    this->ThreadVectorCache = Memory::AllocateArray<Vector3D*>( this->m_NumberOfThreads );
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      this->ThreadVectorCache[thread] = Memory::AllocateArray<Vector3D>( this->ReferenceDims[0] );
  }

  /** Destructor.
   * Free all per-thread data structures.
   */
  virtual ~ParallelElasticFunctional() 
  {
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      if ( ThreadVectorCache[thread] ) 
	Memory::DeleteArray( this->ThreadVectorCache[thread] );
    Memory::DeleteArray( this->ThreadVectorCache );
    
    for ( size_t task = 0; task < this->m_NumberOfThreads; ++task )
      delete this->TaskMetric[task];
    Memory::DeleteArray( this->TaskMetric );
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
	ThreadWarp[thread] = SplineWarpXform::SmartPtr::Null;
	}
      }
  }
  
  /** Evaluate functional after change of a single parameter.
   *@param warp The current deformation.
   *@param metric The metric computed for the base-deformed volume.
   *@param voi Volume-of-Influence for the parameter under consideration.
   *@return The metric after recomputation over the given volume-of-influence.
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
    threadPool.Run( EvaluateGradientThread, InfoTaskGradient );
    
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
      this->WarpedVolume = Memory::AllocateArray<typename VM::Exchange>(  this->DimsX * this->DimsY * this->DimsZ  );

    const size_t numberOfTasks = std::min<size_t>( this->m_NumberOfTasks, this->DimsY * this->DimsZ );
    for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
      {
      InfoTaskComplete[taskIdx].thisObject = this;
      }
    
    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->TaskMetric[taskIdx]->Reset();
      }
    
    ThreadPool::GetGlobalThreadPool().Run( EvaluateCompleteThread, this->InfoTaskComplete );
    
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
};

//@}

} // namespace cmtk

#endif // __cmtkParallelElasticFunctional_h_included_
