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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkParallelElasticFunctional_h_included_
#define __cmtkParallelElasticFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkVoxelMatchingElasticFunctional.h>

#ifdef CMTK_BUILD_SMP
#  include <cmtkThreads.h>
#  include <cmtkMutexLock.h>
#  include <cmtkTimers.h>
#endif

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
template<class VM, class W> 
class ParallelElasticFunctional
  /// Inherit from non-parallel functional.
  : public VoxelMatchingElasticFunctional_Template<VM,W> 
{
protected:
  /// Array of warp transformation objects for the parallel threads.
  SmartPointer<W> *ThreadWarp;

  /// Array of storage for simultaneously retrieving multiple deformed vectors.
  Vector3D **ThreadVectorCache;

  /** Number of threads that this object was created for.
   * All duplicated data structures are generated with the multiplicity given
   * by this value. It is determined from Threads when the object is first
   * instanced. It cannot be changed afterwards.
   */
  int MyNumberOfThreads;

public:
  /// This class.
  typedef ParallelElasticFunctional<VM,W> Self;

  /// Superclass.
  typedef VoxelMatchingElasticFunctional_Template<VM,W> Superclass;

  /// Constructor.
  ParallelElasticFunctional ( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating ) :
    VoxelMatchingElasticFunctional_Template<VM,W>( reference, floating )
  {
    MyNumberOfThreads = Threads::GetNumberOfThreads();

    ThreadWarp = Memory::AllocateArray<typename W::SmartPtr>( MyNumberOfThreads );
    
    InfoThreadGradient = Memory::AllocateArray<typename Self::EvaluateGradientThreadInfo>( MyNumberOfThreads );
    InfoThreadComplete = Memory::AllocateArray<typename Self::EvaluateCompleteThreadInfo>( MyNumberOfThreads );
    
    ThreadMetric = Memory::AllocateArray<VM*>( MyNumberOfThreads );
    for ( int thread = 0; thread < MyNumberOfThreads; ++thread )
      ThreadMetric[thread] = new VM( *(this->Metric) );
    
    ThreadConsistencyHistogram = Memory::AllocateArray<JointHistogram<unsigned int>*>( MyNumberOfThreads );
    for ( int thread = 0; thread < MyNumberOfThreads; ++thread )
      ThreadConsistencyHistogram[thread] = new JointHistogram<unsigned int>;
    
    ThreadVectorCache = Memory::AllocateArray<Vector3D*>( MyNumberOfThreads );
    for ( int thread = 0; thread < MyNumberOfThreads; ++thread )
      ThreadVectorCache[thread] = Memory::AllocateArray<Vector3D>( this->ReferenceDims[0] );
  }

  /** Destructor.
   * Free all per-thread data structures.
   */
  virtual ~ParallelElasticFunctional() 
  {
    for ( int thread = 0; thread < MyNumberOfThreads; ++thread )
      if ( ThreadVectorCache[thread] ) delete[] ThreadVectorCache[thread];
    delete[] ThreadVectorCache;
    
    for ( int thread = 0; thread < MyNumberOfThreads; ++thread )
      delete ThreadMetric[thread];
    delete[] ThreadMetric;
    
    for ( int thread = 0; thread < MyNumberOfThreads; ++thread )
      delete ThreadConsistencyHistogram[thread];
    delete[] ThreadConsistencyHistogram;
    
    delete[] ThreadWarp;
    delete[] InfoThreadGradient;
    delete[] InfoThreadComplete;
  }

  /** Set warp transformation.
   * In the multi-threaded implementation, Warp[0] will be linked directly to
   * the given warp, while for all other threads a copy of the original object
   * is created by a call to WarpXform::Clone().
   */
  virtual void SetWarpXform ( WarpXform::SmartPtr& warp ) 
  {
    this->Superclass::SetWarpXform( warp );
    
    for ( int thread = 0; thread < MyNumberOfThreads; ++thread ) 
      {
      if ( this->Warp ) 
	{
	if ( thread ) 
	  {
	  ThreadWarp[thread] = typename W::SmartPtr( dynamic_cast<W*>( this->Warp->Clone() ) );
	  ThreadWarp[thread]->RegisterVolume( this->ReferenceGrid );
	  } 
	else 
	  {
	  ThreadWarp[thread] = W::SmartPtr::DynamicCastFrom( this->Warp );
	  }
	} 
      else
	{
	ThreadWarp[thread] = W::SmartPtr::Null;
	}
      }
  }
  
  /** Evaluate functional for the complete image data.
   * This function builds the pre-computed deformed floating image that is 
   * later used for rapid gradient computation.
   */
  typename Self::ReturnType EvaluateComplete ( CoordinateVector& v ) 
  {
    this->Metric->Reset();
    if ( ! this->WarpedVolume ) 
      this->WarpedVolume = Memory::AllocateArray<typename VM::Exchange>(  this->DimsX * this->DimsY * this->DimsZ  );
    
    ThreadWarp[0]->SetParamVector( v );

    int numberOfThreads = std::min( MyNumberOfThreads, this->DimsY * this->DimsZ );
    for ( int threadIdx = 0; threadIdx < numberOfThreads; ++threadIdx ) 
      {
      InfoThreadComplete[threadIdx].thisObject = this;
      InfoThreadComplete[threadIdx].ThisThreadIndex = threadIdx;
      InfoThreadComplete[threadIdx].NumberOfThreads = numberOfThreads;
      }
    
    Threads::RunThreads( EvaluateCompleteThread, numberOfThreads, InfoThreadComplete );
    
    return this->WeightedTotal( this->Metric->Get(), ThreadWarp[0] );
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
    
    const typename VM::Exchange unsetY = this->Metric->DataY.padding();
    localMetric->CopyUnsafe( *this->Metric );
    r = voi->startX + this->DimsX * ( voi->startY + this->DimsY * voi->startZ );
    for ( pZ = voi->startZ; pZ<voi->endZ; ++pZ ) 
      {
      for ( pY = voi->startY; pY<voi->endY; ++pY ) 
	{
	pVec = vectorCache;
	warp->GetTransformedGridSequenceNonVirtual( pVec,voi->endX-voi->startX, voi->startX, pY, pZ );
	for ( pX = voi->startX; pX<voi->endX; ++pX, ++r, ++pVec ) 
	  {
	  // Remove this sample from incremental metric according to "ground warp" image.
	  const typename VM::Exchange sampleX = this->Metric->GetSampleX( r );
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
	    else
	      {
	      localMetric->Increment( sampleX, unsetY );
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
    const typename Self::ReturnType current = this->EvaluateComplete( v );

    if ( this->AdaptiveFixParameters && this->WarpNeedsFixUpdate ) 
      {
      this->UpdateWarpFixedParameters();
      }
    
    // Make sure we don't create more threads than we have parameters.
    // Actually, we shouldn't create more than the number of ACTIVE parameters.
    // May add this at some point. Anyway, unless we have A LOT of processors,
    // we shouldn't really ever have more threads than active parameters :))
    int numberOfThreads = std::min<int>( MyNumberOfThreads, this->Dim );

    for ( int threadIdx = 0; threadIdx < numberOfThreads; ++threadIdx ) 
      {
      InfoThreadGradient[threadIdx].thisObject = this;
      InfoThreadGradient[threadIdx].ThisThreadIndex = threadIdx;
      InfoThreadGradient[threadIdx].NumberOfThreads = numberOfThreads;
      InfoThreadGradient[threadIdx].Step = step;
      InfoThreadGradient[threadIdx].Gradient = g.Elements;
      InfoThreadGradient[threadIdx].BaseValue = current;
      InfoThreadGradient[threadIdx].Parameters = &v;
      }

    Threads::RunThreads( EvaluateGradientThread, numberOfThreads, InfoThreadGradient );
    
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

    int numberOfThreads = std::min( MyNumberOfThreads, this->DimsY * this->DimsZ );

    for ( int threadIdx = 0; threadIdx < numberOfThreads; ++threadIdx ) 
      {
      InfoThreadComplete[threadIdx].thisObject = this;
      InfoThreadComplete[threadIdx].ThisThreadIndex = threadIdx;
      InfoThreadComplete[threadIdx].NumberOfThreads = numberOfThreads;
      }
    
    Threads::RunThreads( EvaluateCompleteThread, numberOfThreads, InfoThreadComplete );
    
    return this->WeightedTotal( this->Metric->Get(), ThreadWarp[0] );
  }

private:
  /** Metric object for threadwise computation.
   * The objects in this array are the per-thread equivalent of the
   * VoxelMatchingElasticFunctional::IncrementalMetric object.
   */
  VM** ThreadMetric;
  
  /// Consistency histogram objects for threadwise computation.
  JointHistogram<unsigned int>** ThreadConsistencyHistogram;
  
  /** Thread parameter block for incremental gradient computation.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  class EvaluateGradientThreadInfo 
  {
  public:
    /** Pointer to the functional object that created the thread. */
    Self *thisObject;
    /// Current parameter vector.
    CoordinateVector *Parameters;
    /// Unique index of this thread instance among all threads.
    int ThisThreadIndex;
    /// Total number of threads created.
    int NumberOfThreads;
    /// Current global coordinate stepping.
    typename Self::ParameterType Step;
    /// Pointer to gradient vector that is the target for computation results.
    Types::Coordinate *Gradient;
    /// Base functional value used for comparing new values to.
    double BaseValue;
  };
  
  /// Info blocks for parallel threads evaluating functional gradient.
  typename Self::EvaluateGradientThreadInfo *InfoThreadGradient;
  
  /** Compute functional gradient as a thread.
   * This function (i.e., each thread) iterates over all parameters of the
   * current warp transformation. Among all active (i.e., not disabled)
   * parameters, it selects the ones that have an index with modulus
   * equal to the threads index when divided by the total number of threads.
   * For these parameters, the thread computes the partial derivative of the
   * functional by finite-difference approximation.
   */
  static CMTK_THREAD_RETURN_TYPE EvaluateGradientThread( void* arg ) 
  {
    typename Self::EvaluateGradientThreadInfo *info = static_cast<typename Self::EvaluateGradientThreadInfo*>( arg );
    
    Self *me = info->thisObject;
    SmartPointer<W>& Warp = me->ThreadWarp[info->ThisThreadIndex];
    
    // Set parameter vector on all but first warp (was set on first in 
    // EvaluateComplete() already.
    if ( info->ThisThreadIndex ) Warp->SetParamVector( *info->Parameters );
    
    VM* threadMetric = me->ThreadMetric[info->ThisThreadIndex];
    Vector3D *vectorCache = me->ThreadVectorCache[info->ThisThreadIndex];
    Types::Coordinate *p = Warp->m_Parameters;
    
    Types::Coordinate pOld;
    double upper, lower;
    Rect3D *voi = me->VolumeOfInfluence;

    int activeDims = 0;
    for ( size_t dim = 0; dim < me->Dim; ++dim, ++voi ) 
      {
      if ( me->StepScaleVector[dim] <= 0 ) 
	{
	// let the "0" thread do the common work
	if ( info->ThisThreadIndex == 0 ) info->Gradient[dim] = 0;
	}
      else
	{
	if ( (activeDims % info->NumberOfThreads) == info->ThisThreadIndex ) 
	  {
	  const typename Self::ParameterType thisStep = info->Step * me->StepScaleVector[dim];

	  pOld = p[dim];

	  p[dim] += thisStep;
	  upper = me->EvaluateIncremental( Warp, threadMetric, voi, vectorCache );
	  p[dim] = pOld - thisStep;
	  lower = me->EvaluateIncremental( Warp, threadMetric, voi, vectorCache );
	  
	  p[dim] = pOld;
	  me->WeightedDerivative( lower, upper, Warp, dim, thisStep );
	  
	  if ( (upper > info->BaseValue ) || (lower > info->BaseValue) ) 
	    {
	    // actually, we should divide here by step*StepScaleVector[dim],
	    // shouldn't we?!
	    info->Gradient[dim] = upper - lower;
	    } 
	  else
	    {
	    info->Gradient[dim] = 0;
	    }
	  }
	++activeDims;
	}
      }
    return CMTK_THREAD_RETURN_VALUE;
  }

  /** Thread parameter block for complete functional evaluation.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  class EvaluateCompleteThreadInfo 
  {
  public:
    /** Pointer to the functional object that created the thread. */
    Self *thisObject;
    /// Unique index of this thread instance among all threads.
    int ThisThreadIndex;
    /// Total number of threads created.
    int NumberOfThreads;
  };
  
  /** Info blocks for parallel threads evaluating complete functional. */
  typename Self::EvaluateCompleteThreadInfo *InfoThreadComplete;
    
  /// Multi-threaded implementation of complete metric evaluation.
  static CMTK_THREAD_RETURN_TYPE EvaluateCompleteThread ( void *arg ) 
  {
    typename Self::EvaluateCompleteThreadInfo *info = static_cast<typename Self::EvaluateCompleteThreadInfo*>( arg );
    
    Self *me = info->thisObject;
    const W *warp = me->ThreadWarp[0];
    VM* threadMetric = me->ThreadMetric[info->ThisThreadIndex];
    threadMetric->Reset();
    Vector3D *vectorCache = me->ThreadVectorCache[info->ThisThreadIndex];
    
    typename VM::Exchange* warpedVolume = me->WarpedVolume;
    const typename VM::Exchange unsetY = me->Metric->DataY.padding();
    
    Vector3D *pVec;
    int pX, pY, pZ;
    
    int fltIdx[3];
    Types::Coordinate fltFrac[3];
    
    int rowCount = ( me->DimsY * me->DimsZ );
    int rowFrom = ( rowCount / info->NumberOfThreads ) * info->ThisThreadIndex;
    int rowTo = ( info->ThisThreadIndex == (info->NumberOfThreads-1) ) ? rowCount : ( rowCount / info->NumberOfThreads ) * ( info->ThisThreadIndex + 1 );
    int rowsToDo = rowTo - rowFrom;
    
    int pYfrom = rowFrom % me->DimsY;
    int pZfrom = rowFrom / me->DimsY;
    
    int offset, r = rowFrom * me->DimsX;
    for ( pZ = pZfrom; (pZ < me->DimsZ) && rowsToDo; ++pZ ) 
      {
      for ( pY = pYfrom; (pY < me->DimsY) && rowsToDo; pYfrom = 0, ++pY, --rowsToDo ) 
	{
	warp->GetTransformedGridSequenceNonVirtual( vectorCache, me->DimsX, 0, pY, pZ );
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
	    } 
	  else 
	    {
	    warpedVolume[r] = unsetY;
	    }
	  
	  threadMetric->Increment( me->Metric->GetSampleX(r), warpedVolume[r] );
	  }
	}
      }
    
    me->MetricMutex.Lock();
    me->Metric->AddHistogram( *threadMetric );
    me->MetricMutex.Unlock();
    
    return CMTK_THREAD_RETURN_VALUE;
  }
  
  /// Multi-threaded implementation of complete metric evaluation.
  static CMTK_THREAD_RETURN_TYPE EvaluateThread ( void *const arg ) 
  {
    typename Self::EvaluateCompleteThreadInfo *info = static_cast<typename Self::EvaluateCompleteThreadInfo*>( arg );
    
    Self *const me = info->thisObject;
    W *const warp = me->ThreadWarp[0];
    VM *const threadMetric = me->ThreadMetric[info->ThisThreadIndex];
    VM *const metric = me->Metric;
    threadMetric->Reset();
    Vector3D *const vectorCache = me->ThreadVectorCache[info->ThisThreadIndex];
    
    const typename VM::Exchange unsetY = metric->DataY.padding();
    
    const int dimsX = me->DimsX;
    const int dimsY = me->DimsY;
    const int dimsZ = me->DimsZ;
    
    const int fltDimsX = me->FltDimsX;
    const int fltDimsY = me->FltDimsY;
     
    Vector3D *pVec;
    int pX, pY, pZ;
    int fltIdx[3];
    Types::Coordinate fltFrac[3];
     
    const int rowCount = ( dimsY * dimsZ );
    const int rowFrom = 
      ( rowCount / info->NumberOfThreads ) * info->ThisThreadIndex;
    const int rowTo = ( info->ThisThreadIndex == (info->NumberOfThreads-1) ) 
      ? rowCount
      : ( rowCount / info->NumberOfThreads ) * ( info->ThisThreadIndex + 1 );
    int rowsToDo = rowTo - rowFrom;
     
    const int pYfrom = rowFrom % dimsY;
    const int pZfrom = rowFrom / dimsY;
     
    int offset, r = rowFrom * dimsX;
     
    for ( pZ = pZfrom; (pZ < dimsZ) && rowsToDo; ++pZ ) 
      {
      for ( pY = pYfrom; (pY < dimsY) && rowsToDo; pYfrom = 0, ++pY, --rowsToDo )
	{
	warp->GetTransformedGridSequenceNonVirtual( vectorCache, dimsX, 0, pY, pZ );

	pVec = vectorCache;
	for ( pX = 0; pX < dimsX; ++pX, ++r, ++pVec ) 
	  {
	  // Tell us whether the current location is still within the floating volume and get the respective voxel.
	  Vector3D::CoordMultInPlace( *pVec, me->FloatingInverseDelta );
	  if ( me->FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	    {
	    // Compute data index of the floating voxel in the floating volume.
	    offset = fltIdx[0] + fltDimsX * ( fltIdx[1]+ fltDimsY * fltIdx[2] );
	     
	    // Continue metric computation.
	    threadMetric->Increment( metric->GetSampleX( r ), metric->GetSampleY( offset, fltIdx, fltFrac ) );
	    } 
	  else
	    {
	    threadMetric->Increment( metric->GetSampleX( r ), unsetY );
	    }
	  }
	}
      }
     
    me->MetricMutex.Lock();
    me->Metric->AddHistogram( threadMetric );
    me->MetricMutex.Unlock();
     
    return CMTK_THREAD_RETURN_VALUE;
  }
  
  /** Thread parameter block.
    */
  class UpdateFixedGrayThreadInfo 
  {
  public:
    /// Pointer to the functional object that created the thread.
    Self *thisObject;
    /// Unique index of this thread instance among all threads.
    int ThisThreadIndex;
    /// Total number of threads created.
    int NumberOfThreads;
    /// Maximum entropy of the reference image for this thread.
    double EntropyMaxRef;
    /// Minimum entropy of the reference image for this thread.
    double EntropyMinRef;
    /// Maximum entropy of the floating image for this thread.
    double EntropyMaxFlt;
    /// Minimum entropy of the floating image for this thread.
    double EntropyMinFlt;
  };

  /// Info blocks for parallel threads evaluating complete functional.
  typename Self::UpdateFixedGrayThreadInfo *UpdateFixedGrayThreadInfo;

  static CMTK_THREAD_RETURN_TYPE UpdateWarpFixedParametersGrayThread( void *const arg )
  {
    typename Self::UpdateFixedGrayThreadInfo *info = static_cast<typename Self::UpdateFixedGrayThreadInfo*>( arg );

    Self *me = info->thisObject;
    W *warp = me->ThreadWarp[0];
    VM *threadMetric = me->ThreadMetric[info->ThisThreadIndex];
    VM *metric = me->Metric;

    JointHistogram<unsigned int> *threadHistogram = me->ThreadConsistencyHistogram[info->ThisThreadIndex];
    if ( ! threadHistogram ) 
      {
      threadHistogram = me->ThreadConsistencyHistogram[info->ThisThreadIndex] = new JointHistogram<unsigned int>();

      unsigned int numSamplesX = metric->GetNumberOfSamplesX();
      Types::DataItem fromX, toX;
      metric->GetRangeX( fromX, toX );
      unsigned int numBinsX = threadHistogram->CalcNumBins( numSamplesX, fromX, toX );

      unsigned int numSamplesY = metric->GetNumberOfSamplesY();
      Types::DataItem fromY, toY;
      metric->GetRangeY( fromY, toY );
      unsigned int numBinsY = threadHistogram->CalcNumBins( numSamplesY, fromY, toY );
       
      threadHistogram->SetNumBins( numBinsX, numBinsY );
      threadHistogram->SetRangeX( fromX, toX );
      threadHistogram->SetRangeY( fromY, toY );
      }
     
    const int numCtrlPoints = me->Dim / 3;
     
    Rect3D voi;
    Vector3D fromVOI, toVOI;
    int pX, pY, pZ;

    int inactive = 0;

    double refEntropyMin = HUGE_VAL, refEntropyMax = -HUGE_VAL;
    double fltEntropyMin = HUGE_VAL, fltEntropyMax = -HUGE_VAL;
    double entropyRef, entropyFlt;

    const typename VM::Exchange unsetY = me->Metric->DataY.padding();

    for ( int ctrl = info->ThisThreadIndex; ctrl < numCtrlPoints; ctrl += info->NumberOfThreads ) 
      {
      threadHistogram->Reset();
       
      /// We cannot use the precomputed table of VOIs here because in "fast"
	/// mode, these VOIs are smaller than we want them here.
	warp->GetVolumeOfInfluence( 3 * ctrl, me->ReferenceFrom, me->ReferenceTo, fromVOI, toVOI, 0 );
	me->GetReferenceGridRange( fromVOI, toVOI, voi );
	 
	int r = voi.startX + me->DimsX * ( voi.startY + me->DimsY * voi.startZ );
	 
	const int endLineIncrement = ( voi.startX + ( me->DimsX - voi.endX ) );
	const int endPlaneIncrement = me->DimsX * ( voi.startY + ( me->DimsY - voi.endY ) );
	 
	for ( pZ = voi.startZ; pZ < voi.endZ; ++pZ ) 
	  {
	  for ( pY = voi.startY; pY < voi.endY; ++pY ) 
	    {
	    for ( pX = voi.startX; pX < voi.endX; ++pX, ++r ) 
	      {
	      // Continue metric computation.
	      //	      if ( me->WarpedVolume[r] != unsetY ) {
	      threadHistogram->Increment( metric->GetSampleX( r ), me->WarpedVolume[r] );
	      //	      }
	      }
	    r += endLineIncrement;
	    }
	  r += endPlaneIncrement;
	  }
	threadHistogram->GetMarginalEntropies( entropyRef, entropyFlt );
	 
	if ( entropyRef < refEntropyMin ) refEntropyMin = entropyRef;
	if ( entropyRef > refEntropyMax ) refEntropyMax = entropyRef;
	if ( entropyFlt < fltEntropyMin ) fltEntropyMin = entropyFlt;
	if ( entropyFlt > fltEntropyMax ) fltEntropyMax = entropyFlt;
      }
     
    info->EntropyMinRef = refEntropyMin;
    info->EntropyMaxRef = refEntropyMax;
    info->EntropyMinFlt = fltEntropyMin;
    info->EntropyMaxFlt = fltEntropyMax;
     
    return CMTK_THREAD_RETURN_VALUE;
  }
  
  /** Thread parameter block.
    */
  class UpdateFixedLabelsThreadInfo 
  {
  public:
    /// Pointer to the functional object that created the thread.
    Self *thisObject;
    /// Unique index of this thread instance among all threads.
    int ThisThreadIndex;
    /// Total number of threads created.
    int NumberOfThreads;
  };
  
  /// Info blocks for parallel threads evaluating complete functional.
  typename Self::UpdateFixedLabelsThreadInfo *UpdateFixedLabelsThreadInfo;

  static CMTK_THREAD_RETURN_TYPE UpdateWarpFixedParametersLabelsThread( void *const arg )
  {
    typename Self::UpdateFixedLabelsThreadInfo *info = static_cast<typename Self::UpdateFixedLabelsThreadInfo*>( arg );
     
    Self *me = info->thisObject;
    W *warp = me->ThreadWarp[0];
    VM *metric = me->Metric;

    const int numCtrlPoints = me->Dim / 3;

    Rect3D voi;
    Vector3D fromVOI, toVOI;
    int pX, pY, pZ;
     
    int inactive = 0;
     
    for ( int ctrl = info->ThisThreadIndex; ctrl < numCtrlPoints; ctrl += info->NumberOfThreads ) 
      {
      /// We cannot use the precomputed table of VOIs here because in "fast"
       /// mode, these VOIs are smaller than we want them here.
      warp->GetVolumeOfInfluence( 3 * ctrl, me->ReferenceFrom, me->ReferenceTo, fromVOI, toVOI, 0 );
      me->GetReferenceGridRange( fromVOI, toVOI, voi );
       
      int r = voi.startX + me->DimsX * ( voi.startY + me->DimsY * voi.startZ );
       
      const int endLineIncrement = ( voi.startX + (me->DimsX-voi.endX) );
      const int endPlaneIncrement = me->DimsX * ( voi.startY + (me->DimsY-voi.endY) );
       
      bool active = false;
      for ( pZ = voi.startZ; (pZ < voi.endZ) && !active; ++pZ ) 
	{
	for ( pY = voi.startY; (pY < voi.endY) && !active; ++pY ) 
	  {
	  for ( pX = voi.startX; (pX < voi.endX); ++pX, ++r ) 
	    {
	    // Continue metric computation.
	    if ( ( metric->GetSampleX( r ) != 0 ) || ( ( me->WarpedVolume[r] != VM::unset() ) && ( me->WarpedVolume[r] != 0 ) ) ) 
	      {
	      active = true;
	      break;
	      }
	    }
	  r += endLineIncrement;
	  }
	r += endPlaneIncrement;
	}
       
      if ( !active ) 
	{
	inactive += 3;
	 
	int dim = 3 * ctrl;
	for ( int idx=0; idx<3; ++idx, ++dim ) 
	  {
	  warp->SetParameterInactive( dim );
	  me->StepScaleVector[dim] = me->GetParamStep( dim );
	  }
	}
      }
     
    return CMTK_THREAD_RETURN_VALUE;
  }
};

//@}

} // namespace cmtk

#endif // __cmtkParallelElasticFunctional_h_included_
