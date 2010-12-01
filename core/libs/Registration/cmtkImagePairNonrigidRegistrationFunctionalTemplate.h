/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkImagePairNonrigidRegistrationFunctionalTemplate_h_included_
#define __cmtkImagePairNonrigidRegistrationFunctionalTemplate_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkImagePairNonrigidRegistrationFunctional.h>

#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkDataTypeTraits.h>

#ifdef CMTK_BUILD_DEMO
#  include <IO/cmtkXformIO.h>
#endif // #ifdef CMTK_BUILD_DEMO

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
template<class VM> 
class ImagePairNonrigidRegistrationFunctionalTemplate
  /// Inherit from general image pair registration functional.
  : public ImagePairNonrigidRegistrationFunctional 
{
protected:
  /** Metric object for incremental computation.
   * Before computing the incremental metric after change of one parameter,
   * the global metric is copied to this object. It is then used for in-place
   * application of all necessary changes, leaving the original metric intact.
   *@see #EvaluateIncremental
   */
  SmartPointer<VM> m_IncrementalMetric;
  
public:
  /// This class.
  typedef ImagePairNonrigidRegistrationFunctionalTemplate<VM> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef ImagePairNonrigidRegistrationFunctional Superclass;

  /// Constructor.
  ImagePairNonrigidRegistrationFunctionalTemplate<VM>( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating, const Interpolators::InterpolationEnum interpolation )
    : ImagePairNonrigidRegistrationFunctional( reference, floating )
  {
    this->m_InfoTaskGradient.resize( this->m_NumberOfTasks );
    this->m_InfoTaskComplete.resize( this->m_NumberOfTasks );
    
    this->m_Metric = ImagePairSimilarityMeasure::SmartPtr( new VM( reference, floating, interpolation ) );
    this->m_TaskMetric.resize( this->m_NumberOfThreads, dynamic_cast<const VM&>( *this->m_Metric ) );
  }

  /** Destructor.
   * Free all per-thread data structures.
   */
  virtual ~ImagePairNonrigidRegistrationFunctionalTemplate<VM>() {}

  /// Set flag and value for forcing values outside the floating image.
  virtual void SetForceOutside
  ( const bool flag = true, const Types::DataItem value = 0 )
  {
    this->m_ForceOutsideFlag = flag;
    this->m_ForceOutsideValueRescaled = this->m_Metric->GetFloatingValueScaled( value );
  }

  /** Set warp transformation.
   */
  virtual void SetWarpXform ( SplineWarpXform::SmartPtr& warp )
  {
    Superclass::SetWarpXform( warp );
    this->WarpNeedsFixUpdate = true;
  }
  
  /// Match intensities of reference and floating images.
  virtual void MatchRefFltIntensities();

  /** Evaluate functional for the complete image data.
   * This function builds the pre-computed deformed floating image that is 
   * later used for rapid gradient computation.
   */
  typename Self::ReturnType Evaluate() 
  {
    this->m_Metric->Reset();
    if ( ! this->m_WarpedVolume )
      {
      this->m_WarpedVolume = Memory::AllocateArray<Types::DataItem>(  this->m_DimsX * this->m_DimsY * this->m_DimsZ  );
      }
    
    const size_t numberOfTasks = std::min<size_t>( this->m_NumberOfTasks, this->m_DimsY * this->m_DimsZ );
    for ( size_t taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
      {
      this->m_InfoTaskComplete[taskIdx].thisObject = this;
      }

    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      this->m_TaskMetric[taskIdx].Reset();
      }
    
    ThreadPool::GetGlobalThreadPool().Run( EvaluateCompleteThread, this->m_InfoTaskComplete );
    
    for ( size_t taskIdx = 0; taskIdx < this->m_NumberOfThreads; ++taskIdx ) 
      {
      dynamic_cast<VM&>( *(this->m_Metric) ).Add( this->m_TaskMetric[taskIdx] );
      }
    
    return this->WeightedTotal( this->m_Metric->Get(), *(this->m_ThreadWarp[0]) );
  }

  /** Evaluate functional after change of a single parameter.
   *@param warp The current deformation.
   *@param metric The metric computed for the base-deformed volume.
   *@param voi Volume-of-Influence for the parameter under consideration.
   *@return The metric after recomputation over the given volume-of-influence.
   */
  typename Self::ReturnType EvaluateIncremental( const SplineWarpXform& warp, VM& localMetric, const DataGrid::RegionType& voi, Vector3D *const vectorCache ) 
  {
    Vector3D *pVec;
    int pX, pY, pZ, r;
    int fltIdx[3];
    Types::Coordinate fltFrac[3];

    int endLineIncrement = ( voi.From()[0] + ( this->m_DimsX - voi.To()[0]) );
    int endPlaneIncrement = this->m_DimsX * ( voi.From()[1] + (this->m_DimsY - voi.To()[1]) );
    
    const Types::DataItem unsetY = DataTypeTraits<Types::DataItem>::ChoosePaddingValue();
    localMetric = dynamic_cast<VM&>( *this->m_Metric );
    r = voi.From()[0] + this->m_DimsX * ( voi.From()[1] + this->m_DimsY * voi.From()[2] );
    for ( pZ = voi.From()[2]; pZ<voi.To()[2]; ++pZ ) 
      {
      for ( pY = voi.From()[1]; pY<voi.To()[1]; ++pY ) 
	{
	pVec = vectorCache;
	warp.GetTransformedGridRow( voi.To()[0]-voi.From()[0], pVec, voi.From()[0], pY, pZ );
	for ( pX = voi.From()[0]; pX<voi.To()[0]; ++pX, ++r, ++pVec ) 
	  {
	  // Remove this sample from incremental metric according to "ground warp" image.
	  Types::DataItem sampleX;
	  if ( this->m_Metric->GetSampleX( sampleX, r ) )
	    {
	    if ( this->m_WarpedVolume[r] != unsetY )
	      localMetric.Decrement( sampleX, this->m_WarpedVolume[r] );
	    
	    // Tell us whether the current location is still within the floating volume and get the respective voxel.
	    *pVec *= this->m_FloatingInverseDelta;
	    if ( this->m_FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	      {
	      // Continue metric computation.
	      localMetric.Increment( sampleX, this->m_Metric->GetSampleY( fltIdx, fltFrac ) );
	      } 
	    else
	      {
	      if ( this->m_ForceOutsideFlag )
		{
		localMetric.Increment( sampleX, this->m_ForceOutsideValueRescaled );
		}
	      }
	    }
	  }
	r += endLineIncrement;
	}
      r += endPlaneIncrement;
      }
    
    return localMetric.Get();
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
      this->m_InfoTaskGradient[taskIdx].thisObject = this;
      this->m_InfoTaskGradient[taskIdx].Step = step;
      this->m_InfoTaskGradient[taskIdx].Gradient = g.Elements;
      this->m_InfoTaskGradient[taskIdx].BaseValue = current;
      this->m_InfoTaskGradient[taskIdx].Parameters = &v;
      }

    ThreadPool::GetGlobalThreadPool().Run( EvaluateGradientThread, this->m_InfoTaskGradient );
    
    return current;
  }

  /// Evaluate functional.
  virtual typename Self::ReturnType EvaluateAt ( CoordinateVector& v )
  {
    this->m_ThreadWarp[0]->SetParamVector( v );
    return this->Evaluate();
  }

#ifdef CMTK_BUILD_DEMO
  /// Create a snapshot (to disk) of current functional result.
  virtual void SnapshotAt( ParameterVectorType& v )
  {
    this->m_ThreadWarp[0]->SetParamVector( v );
    static int it = 0;
    char path[PATH_MAX];
    snprintf( path, PATH_MAX, "warp-%03d.xform", it++ );
    XformIO::Write( this->m_ThreadWarp[0], path );
  }
#endif

private:
  /** Metric object for threadwise computation.
   * The objects in this array are the per-thread equivalent of the
   * ImagePairNonrigidRegistrationFunctional::IncrementalMetric object.
   */
  std::vector<VM> m_TaskMetric;
  
  /// Consistency histogram objects for threadwise computation.
  JointHistogram<unsigned int>** m_ThreadConsistencyHistogram;
  
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
  std::vector<typename Self::EvaluateGradientTaskInfo> m_InfoTaskGradient;
  
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

    SplineWarpXform& myWarp = *(me->m_ThreadWarp[threadIdx]);
    myWarp.SetParamVector( *info->Parameters );
    
    VM& threadMetric = me->m_TaskMetric[threadIdx];
    Vector3D *vectorCache = me->m_ThreadVectorCache[threadIdx];
    Types::Coordinate *p = myWarp.m_Parameters;
    
    Types::Coordinate pOld;
    double upper, lower;

    const DataGrid::RegionType *voi = me->VolumeOfInfluence + taskIdx;
    for ( size_t dim = taskIdx; dim < me->Dim; dim+=taskCnt, voi+=taskCnt ) 
      {
      if ( me->m_StepScaleVector[dim] <= 0 ) 
	{
	info->Gradient[dim] = 0;
	}
      else
	{
	const typename Self::ParameterType thisStep = info->Step * me->m_StepScaleVector[dim];
	
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
  std::vector<typename Self::EvaluateCompleteTaskInfo> m_InfoTaskComplete;
    
  /// Multi-threaded implementation of complete metric evaluation.
  static void EvaluateCompleteThread ( void *arg, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t ) 
  {
    typename Self::EvaluateCompleteTaskInfo *info = static_cast<typename Self::EvaluateCompleteTaskInfo*>( arg );
    
    Self *me = info->thisObject;
    const SplineWarpXform& warp = *(me->m_ThreadWarp[0]);
    VM& threadMetric = me->m_TaskMetric[threadIdx];
    Vector3D *vectorCache = me->m_ThreadVectorCache[threadIdx];
    
    Types::DataItem* warpedVolume = me->m_WarpedVolume;
    const Types::DataItem unsetY = ( me->m_ForceOutsideFlag ) ? me->m_ForceOutsideValueRescaled : DataTypeTraits<Types::DataItem>::ChoosePaddingValue();
    
    Vector3D *pVec;
    int pX, pY, pZ;
    
    int fltIdx[3];
    Types::Coordinate fltFrac[3];
    
    int rowCount = ( me->m_DimsY * me->m_DimsZ );
    int rowFrom = ( rowCount / taskCnt ) * taskIdx;
    int rowTo = ( taskIdx == (taskCnt-1) ) ? rowCount : ( rowCount / taskCnt ) * ( taskIdx + 1 );
    int rowsToDo = rowTo - rowFrom;
    
    int pYfrom = rowFrom % me->m_DimsY;
    int pZfrom = rowFrom / me->m_DimsY;
    
    int r = rowFrom * me->m_DimsX;
    for ( pZ = pZfrom; (pZ < me->m_DimsZ) && rowsToDo; ++pZ ) 
      {
      for ( pY = pYfrom; (pY < me->m_DimsY) && rowsToDo; pYfrom = 0, ++pY, --rowsToDo ) 
	{
	warp.GetTransformedGridRow( me->m_DimsX, vectorCache, 0, pY, pZ );
	pVec = vectorCache;
	for ( pX = 0; pX<me->m_DimsX; ++pX, ++r, ++pVec ) 
	  {
	  // Tell us whether the current location is still within the 
	  // floating volume and get the respective voxel.
	  *pVec *= me->m_FloatingInverseDelta;
	  if ( me->m_FloatingGrid->FindVoxelByIndex( *pVec, fltIdx, fltFrac ) ) 
	    {
	    // Continue metric computation.
	    warpedVolume[r] = me->m_Metric->GetSampleY( fltIdx, fltFrac );
	    
	    Types::DataItem value;
	    if ( me->m_Metric->GetSampleX( value, r ) )
	      {
	      threadMetric.Increment( value, warpedVolume[r] );
	      }
	    } 
	  else 
	    {
	    warpedVolume[r] = unsetY;
	    }
	  }
	}
      }
  }

private:
  /** Warp's fixed parameters need to be updated.
   * This flag is set when the warp transformation is set or modified. It
   * signals that the active and passive parameters of the transformation
   * will have to be updated before the next gradient computation.
   */
  bool WarpNeedsFixUpdate;

  /// Histogram used for consistency computation.
  JointHistogram<unsigned int>::SmartPtr m_ConsistencyHistogram;

  /** Update set of active and passive parameters.
   * This function computes local entropies in the neighborhood of all control
   * points of the Warp transformation. Those control points for which both
   * reference and floating image have less than half the maximum entropy in
   * this neighborhood as compared to the rest of the image are set passive.
   * The passive parameters are not considered for gradient computation and
   * therefore save significant computation time.
   */
  void UpdateWarpFixedParameters();
};

//@}

} // namespace cmtk

#include "cmtkImagePairNonrigidRegistrationFunctionalTemplate.txx"

#endif // __cmtkImagePairNonrigidRegistrationFunctionalTemplate_h_included_
