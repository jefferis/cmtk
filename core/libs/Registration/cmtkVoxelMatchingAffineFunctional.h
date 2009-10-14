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

#ifndef __cmtkVoxelMatchingAffineFunctional_h_included_
#define __cmtkVoxelMatchingAffineFunctional_h_included_

#include <cmtkconfig.h>

#include <assert.h>

#include <cmtkVoxelMatchingFunctional.h>

#include <cmtkVector.h>
#include <cmtkAffineXform.h>
#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>

#include <cmtkMathUtil.h>
#include <cmtkTypes.h>
#include <cmtkException.h>

#include <cmtkVolumeClipping.h>
#include <cmtkVolumeAxesHash.h>

#include <cmtkThreadPool.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base-class for affine registration functionals.
 */
class VoxelMatchingAffineFunctional : 
  /// Inherit from voxel matching functional.
  public VoxelMatchingFunctional 
{
public:
  /// This class type.
  typedef VoxelMatchingAffineFunctional Self;

  /// Superclass.
  typedef VoxelMatchingFunctional Superclass;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v )  
  {
    this->m_AffineXform->GetParamVector( v );
  }

  /// Return parameter stepping.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    return this->m_AffineXform->GetParamStep( idx, FloatingSize, mmStep );
  }

  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const 
  {
    return this->m_AffineXform->ParamVectorDim();
  }

  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const 
  {
    return this->m_AffineXform->VariableParamVectorDim();
  }

protected:
  /// Current coordinate transformation.
  AffineXform::SmartPtr m_AffineXform;

  /// Utility object for volume clipping.
  VolumeClipping Clipper;

  /** Perform clipping/cropping in z-direction.
   * This function computes the intersection of reference and floating data in
   * z-direction. It determines the range of indices of those planes in the
   * reference that intersect the floating. This is the range over which to 
   * for-loop during metric computation.
   *@param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *@param origin Starting point of the reference volume.
   *@param start Upon return, this reference is set to the index of first plane
   * in the reference that intersects the floating.
   *@param end Upon return, this reference is set to one plus the index of the
   * last plane in the reference that intersects the floating.
   *@return 1 if there is an intersection of reference and floating, 0 if there
   * isn't. The range of indices returned in "start" and "end" is only
   * guaranteed to be valid if 1 is the return value.
   */
  int ClipZ ( const VolumeClipping& clipper, const Vector3D& origin, GridIndexType& start, GridIndexType &end ) const
  {
    // perform clipping
    Types::Coordinate fromFactor, toFactor;
    if (! clipper.ClipZ( fromFactor, toFactor, origin ) )
      return 0;

    // there is an intersection: Look up the corresponding grid indices
    start = static_cast<GridIndexType>( (ReferenceDims[2]-1)*fromFactor );
    end = 1+std::min( (int)(ReferenceDims[2]-1), (int)(1 + ((ReferenceDims[2]-1)*toFactor) ) );
    
    // finally, apply cropping boundaries of the reference volume
    start = std::max<GridIndexType>( start, ReferenceCropFrom[2] );
    end = std::min<GridIndexType>( end, ReferenceCropTo[2] );
    
    // return 1 iff index range is non-empty.
    return (start < end );
  }

  /** Perform clipping/cropping in x-direction.
   * This function computes the intersection of reference and floating data in
   * x-direction. It determines the range of indices of those voxels in the
   * current reference row that intersect the floating image. This is the range
   * over which to for-loop during metric computation.
   *
   * Compared to ClipZ and ClipY, this step has to operate very exact as there
   * is no further level that would reduce remaining invalid voxels. Therefore,
   * clipper.ClipX() is called with an extended initial range of indices and an
   * explicitly open upper bound.
   *
   * This is necessary to discriminate inside-boundary from on-boundary voxels.
   * For the right, upper and back boundary, on-boundary voxels are already
   * outside the allowed range as the upper boundaries of the volume are open
   * in terms of interpolation.
   *@param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *@param origin Starting point of the current row in the reference volume.
   *@param start Upon return, this reference is set to the index of first voxel
   * in the reference that intersects the floating image.
   *@param end Upon return, this reference is set to one plus the index of the
   * last voxel in the reference that intersects the floating image.
   *@return 1 if there is an intersection of the current reference row and
   * the floating, 0 if there isn't. The range of indices returned in "start"
   * and "end" is only guaranteed to be valid if 1 is the return value.
   */
  int ClipX ( const VolumeClipping& clipper, const Vector3D& origin, GridIndexType& start, GridIndexType &end ) const
  {
    // perform clipping
    Types::Coordinate fromFactor, toFactor;
    if ( ! clipper.ClipX( fromFactor, toFactor, origin, 0, 2, false, true ) )
      return 0;

    fromFactor = std::min<Types::Coordinate>( 1.0, fromFactor );
	      
    // there is an intersection: Look up the corresponding grid indices
    start = std::max( 0, (int)((ReferenceDims[0]-1)*fromFactor)-1 );
    while ( ( start*ReferenceGrid->m_Delta[0] < fromFactor*ReferenceSize[0]) && ( start < ReferenceDims[0] ) ) 
      ++start;
    
    if ( (toFactor > 1.0) || (start == ReferenceDims[0]) ) 
      {
      end = ReferenceDims[0];
      } 
    else
      {
      end = std::min( ReferenceDims[0]-2, (int)(1 + (ReferenceDims[0]-1)*toFactor));
      while ( end*ReferenceGrid->m_Delta[0] > toFactor*ReferenceSize[0] ) // 'if' not sufficient!	
	--end;
      ++end; // otherwise end=1+min(...) and ...[0][end-1] above!!
      }
    
    // finally, apply cropping boundaries of the reference volume
    start = std::max<GridIndexType>( start, ReferenceCropFrom[0] );
    end = std::min<GridIndexType>( end, ReferenceCropTo[0] );
    
    // return 1 iff index range is non-empty.
    return (start < end );
  }

  /** Perform clipping/cropping in y-direction.
   * This function computes the intersection of reference and floating data in
   * y-direction. It determines the range of indices of those rows in the
   * current reference plane that intersect the floating image. This is the
   * range over which to for-loop during metric computation.
   *@param clipper A volume clipping object with clipping boundaries and grid
   * orientation set.
   *@param origin Starting point of the current plane in the reference volume.
   *@param start Upon return, this reference is set to the index of first row
   * in the reference that intersects the floating image.
   *@param end Upon return, this reference is set to one plus the index of the
   * last row in the reference that intersects the floating image.
   *@return 1 if there is an intersection of the current reference plane and
   * the floating, 0 if there isn't. The range of indices returned in "start" 
   * and "end" is only guaranteed to be valid if 1 is the return value.
   */
  int ClipY ( const VolumeClipping& clipper, const Vector3D& origin, GridIndexType& start, GridIndexType &end ) const
  {
    // perform clipping
    Types::Coordinate fromFactor, toFactor;
    if ( !clipper.ClipY( fromFactor, toFactor, origin ) )
      return 0;

    // there is an intersection: Look up the corresponding grid indices
    start = static_cast<GridIndexType>( (ReferenceDims[1]-1)*fromFactor );
    
    if ( toFactor > 1.0 ) 
      {
      end = ReferenceDims[1];
      } 
    else
      {
      end = 1+std::min( ReferenceDims[1]-1, (int)(1+(ReferenceDims[1]-1)*toFactor ) );
      }
    // finally, apply cropping boundaries of the reference volume
    start = std::max<GridIndexType>( start, ReferenceCropFrom[1] );
    end = std::min<GridIndexType>( end, ReferenceCropTo[1] );
    
    // return 1 iff index range is non-empty.
    return (start < end );
  }

public:
  /// Constructor.
  VoxelMatchingAffineFunctional( UniformVolume::SmartPtr refVolume, UniformVolume::SmartPtr modVolume, AffineXform::SmartPtr& affineXform ) 
    : VoxelMatchingFunctional( refVolume, modVolume ) 
  {
    if ( affineXform.IsNull() ) throw Exception();
    this->m_AffineXform = affineXform;
  }

  /// Destructor.
  virtual ~VoxelMatchingAffineFunctional() {}
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
class VoxelMatchingAffineFunctional_Template :
  /// Inherit from affine voxel matching functional
  public VoxelMatchingAffineFunctional, 
  /// Inherit from metric template functional.
  public VoxelMatchingFunctional_Template<VM> 
{
public:
  /// This class type.
  typedef VoxelMatchingAffineFunctional_Template<VM> Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef VoxelMatchingAffineFunctional Superclass;

  /// Return type.
  typedef Functional::ReturnType ReturnType;

  /** Constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   *@param reference The reference (i.e. static) volume.
   *@param floating The floating (i.e. transformed) volume.
   *@param xform A transformation template. This object determines the type
   * of transformation to be optimized. Its initial value is not relevant.
   */
  VoxelMatchingAffineFunctional_Template( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating, AffineXform::SmartPtr& affineXform )
    : VoxelMatchingAffineFunctional( reference, floating, affineXform ),
      VoxelMatchingFunctional_Template<VM>( reference, floating ) 
  {
    this->m_NumberOfThreads = this->m_ThreadPool.GetNumberOfThreads();
    this->m_NumberOfTasks = 2 * this->m_NumberOfThreads;

    this->m_EvaluateTaskInfo = Memory::AllocateArray<typename Self::EvaluateTaskInfo>( this->m_NumberOfTasks );
    for ( size_t threadIdx = 0; threadIdx < this->m_NumberOfTasks; ++threadIdx ) 
      {
      this->m_EvaluateTaskInfo[threadIdx].thisObject = this;
      }

    this->m_ThreadMetric = Memory::AllocateArray<VM*>( m_NumberOfThreads );
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
      this->m_ThreadMetric[thread] = new VM( *(this->Metric) );
  }

  /// Destructor.
  virtual ~VoxelMatchingAffineFunctional_Template() 
  {
    for ( size_t thread = 0; thread < m_NumberOfThreads; ++thread )
      delete m_ThreadMetric[thread];
    Memory::DeleteArray( this->m_ThreadMetric );
    
    Memory::DeleteArray( this->m_EvaluateTaskInfo );
  }

  /// Evaluate with new parameter vector.
  virtual typename Self::ReturnType EvaluateAt ( CoordinateVector& v ) 
  {
    this->m_AffineXform->SetParamVector( v );
    return this->Evaluate();
  }

  /** Compute functional value with volume clipping.
   * This function iterates over all voxels of the reference image that - after
   * applying the current coordinate transformation - are located inside the
   * mode image. This set of voxels is determined on-the-fly by an extension of
   * Liang and Barsky's "Parameterized Line-Clipping" technique.
   *
   * From the resulting sequence of reference/floating voxel pairs, the 
   * selected voxel-based similarity measure (metric) is computed.
   *@param v The current parameter vector describing the effective coordinate
   * transformation.
   *@return The computed similarity measure as returned by the "Metric" 
   * subobject.
   *@see VolumeClipping
   */
  virtual typename Self::ReturnType Evaluate() 
  {
    const VolumeAxesHash axesHash( *this->ReferenceGrid, this->m_AffineXform, this->FloatingGrid->m_Delta, this->FloatingGrid->m_Origin.XYZ );
    const Vector3D *axesHashX = axesHash[0], *axesHashY = axesHash[1], *axesHashZ = axesHash[2];
    
    this->Metric->Reset();

    const int *Dims = this->ReferenceGrid->GetDims();
    const int DimsX = Dims[0], DimsY = Dims[1], DimsZ = Dims[2];

    this->Clipper.SetDeltaX( axesHashX[DimsX-1] - axesHashX[0] );
    this->Clipper.SetDeltaY( axesHashY[DimsY-1] - axesHashY[0] );
    this->Clipper.SetDeltaZ( axesHashZ[DimsZ-1] - axesHashZ[0] );
    this->Clipper.SetClippingBoundaries( this->FloatingCropFromIndex, this->FloatingCropToIndex );
    
    GridIndexType startZ, endZ;
    if ( this->ClipZ( this->Clipper, axesHashZ[0], startZ, endZ ) ) 
      {
      startZ = std::max<GridIndexType>( startZ, this->ReferenceCropFrom[2] );
      endZ = std::min<GridIndexType>( endZ, this->ReferenceCropTo[2] + 1 );
      
      const int numberOfTasks = std::min<size_t>( this->m_NumberOfTasks, endZ - startZ + 1 );
      
      for ( int taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
	{
	this->m_EvaluateTaskInfo[taskIdx].AxesHash = &axesHash;
	this->m_EvaluateTaskInfo[taskIdx].StartZ = startZ;
	this->m_EvaluateTaskInfo[taskIdx].EndZ = endZ;
	}
      
      for ( size_t threadIdx = 0; threadIdx < this->m_NumberOfThreads; ++threadIdx ) 
	{
	this->m_ThreadMetric[threadIdx]->Reset();
	}

      this->m_ThreadPool.Run( EvaluateThread, numberOfTasks, this->m_EvaluateTaskInfo );

      for ( size_t threadIdx = 0; threadIdx < this->m_NumberOfThreads; ++threadIdx ) 
	{
	this->Metric->AddHistogram( *this->m_ThreadMetric[threadIdx] );
	}
      }
    return this->Metric->Get();
  }
  
  /** Number of threads that this object was created for.
   * This is the actual maximum number of threads running at any time, but not
   * necessarily the number of parallel tasks to be completed.
   * All duplicated data structures are generated with the multiplicity given
   * by this value. It is determined from Threads when the object is first
   * instanced. It cannot be changed afterwards.
   */
  size_t m_NumberOfThreads;

  /// Pool of persistent threads to repeatedly evaluate this functional.
  ThreadPool m_ThreadPool;

  /** Number of tasks executed by the thread pool.
   * This is the total number of tasks that the complete evaluation is broken into.
   * Ideally, this should be larger than the number of threads to allow load balancing.
   */
  size_t m_NumberOfTasks;

  /// Metric objects for the separate threads.
  VM** m_ThreadMetric;

  /** Thread parameter block for incremental gradient computation.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  typedef struct 
  {
    /// Pointer to the functional object that created the thread.
    Self *thisObject;
    /// Unique index of this thread instance among all threads.
    int ThisThreadIndex;
    /// Total number of threads created.
    int NumberOfThreads;
    /// Axes hash.
    const VolumeAxesHash* AxesHash;
    /// First plane of clipped reference volume.
    GridIndexType StartZ;
    /// Last plane of clipped reference volume.
    GridIndexType EndZ;
  } EvaluateTaskInfo;
 
  /// Info blocks for parallel threads evaluating functional gradient.
  typename Self::EvaluateTaskInfo *m_EvaluateTaskInfo;

  /** Compute functional gradient as a thread.
    * This function (i.e., each thread) iterates over all parameters of the
    * current warp transformation. Among all active (i.e., not disabled)
    * parameters, it selects the ones that have an index with modulus
    * equal to the threads index when divided by the total number of threads.
    * For these parameters, the thread computes the partial derivative of the
    * functional by finite-difference approximation.
    */
  static void EvaluateThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCont ) 
  {
  typename Self::EvaluateTaskInfo *info = static_cast<typename Self::EvaluateTaskInfo*>( args );

    Self *me = info->thisObject;
#ifndef CMTK_PVI_HISTOGRAMS
    const VM* Metric = me->Metric;
#endif
    VM* threadMetric = me->m_ThreadMetric[threadIdx];

    const Vector3D *hashX = (*info->AxesHash)[0], *hashY = (*info->AxesHash)[1], *hashZ = (*info->AxesHash)[2];
    Vector3D pFloating;

    const int *Dims = me->ReferenceGrid->GetDims();
    const int DimsX = Dims[0], DimsY = Dims[1];

    int fltIdx[3];
    Types::Coordinate fltFrac[3];

    const int FltDimsX = me->FloatingDims[0], FltDimsY = me->FloatingDims[1];

    Vector3D rowStart;
    Vector3D planeStart;

    int offset;
    GridIndexType pX, pY, pZ;
    // Loop over all remaining planes
    for ( pZ = info->StartZ + taskIdx; pZ < info->EndZ; pZ += taskCnt ) 
      {
      // Offset of current reference voxel
      int r = pZ * DimsX * DimsY;
      
      planeStart = hashZ[pZ];
      
      GridIndexType startY, endY;
      if ( me->ClipY( me->Clipper, planeStart, startY, endY ) ) 
	{	
	startY = std::max<GridIndexType>( startY, me->ReferenceCropFrom[1] );
	endY = std::min<GridIndexType>( endY, me->ReferenceCropTo[1] + 1 );
	r += startY * DimsX;
	
	// Loop over all remaining rows
	for ( pY = startY; pY<endY; ++pY ) 
	  {
	  (rowStart = planeStart) += hashY[pY];
	  
	  GridIndexType startX, endX;
	  if ( me->ClipX( me->Clipper, rowStart, startX, endX ) ) 
	    {	    
	    startX = std::max<GridIndexType>( startX, me->ReferenceCropFrom[0] );
	    endX = std::min<GridIndexType>( endX, me->ReferenceCropTo[0] + 1 );
	    
	    r += startX;
	    // Loop over all remaining voxels in current row
	    for ( pX = startX; pX<endX; ++pX, ++r ) 
	      {
	      (pFloating = rowStart) += hashX[pX];
	      
	      // probe volume and get the respective voxel
	      if ( me->FloatingGrid->FindVoxelByIndex( pFloating, fltIdx, fltFrac ) )
		{
		// Compute data index of the floating voxel in the floating 
		// volume.
		offset = fltIdx[0]+FltDimsX*(fltIdx[1]+FltDimsY*fltIdx[2]);
		
		// Continue metric computation.
#ifdef CMTK_PVI_HISTOGRAMS
		threadMetric->Proceed( r, offset, fltIdx, fltFrac );
#else
		threadMetric->Increment( Metric->GetSampleX( r ), Metric->GetSampleY( offset, fltFrac ) );
#endif
		}
	      }
	    r += (DimsX-endX);
	    } 
	  else
	    {
	    r += DimsX;
	    }
	  }
	
	r += (DimsY-endY) * DimsX;
	} 
      else
	{
	r += DimsY * DimsX;
	}
      }
  }
};

/** Constructor function for affine voxel registration functionals.
 * This function takes the index of a metric in the list of available voxel
 * similarity measures plus all required objects. It the creates an appropriate
 * instance of VoxelMatchingAffineFunctional with the correct metric class as template
 * parameter.
 */
VoxelMatchingAffineFunctional* 
CreateAffineFunctional
( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& modVolume, AffineXform::SmartPtr& affineXform );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVoxelMatchingAffineFunctional_h_included_
