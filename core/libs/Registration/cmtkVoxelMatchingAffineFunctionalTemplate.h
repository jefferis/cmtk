/*
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

#ifndef __cmtkVoxelMatchingAffineFunctionalTemplate_h_included_
#define __cmtkVoxelMatchingAffineFunctionalTemplate_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkVoxelMatchingAffineFunctional.h>
#include <Registration/cmtkVoxelMatchingFunctional.h>

#include <Base/cmtkTransformedVolumeAxes.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

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
class VoxelMatchingAffineFunctionalTemplate :
  /// Inherit from affine voxel matching functional
  public VoxelMatchingAffineFunctional, 
  /// Inherit from metric template functional.
  public VoxelMatchingFunctional_Template<VM> 
{
public:
  /// This class type.
  typedef VoxelMatchingAffineFunctionalTemplate<VM> Self;

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
  VoxelMatchingAffineFunctionalTemplate( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating, AffineXform::SmartPtr& affineXform )
    : VoxelMatchingAffineFunctional( reference, floating, affineXform ),
      VoxelMatchingFunctional_Template<VM>( reference, floating ),
      m_NumberOfThreads( ThreadPool::GetGlobalThreadPool().GetNumberOfThreads() )    
  {
    this->m_ThreadMetric.resize( m_NumberOfThreads, dynamic_cast<const VM&>( *(this->Metric) ) );
  }

  /// Destructor.
  virtual ~VoxelMatchingAffineFunctionalTemplate() {}

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
   *@return The computed similarity measure as returned by the "Metric" 
   * subobject.
   *@see VolumeClipping
   */
  virtual typename Self::ReturnType Evaluate() 
  {
    const TransformedVolumeAxes axesHash( *this->ReferenceGrid, this->m_AffineXform, this->FloatingGrid->Deltas().begin(), this->FloatingGrid->m_Offset.begin() );
    const Vector3D *axesHashX = axesHash[0], *axesHashY = axesHash[1], *axesHashZ = axesHash[2];
    
    this->Metric->Reset();

    const DataGrid::IndexType& Dims = this->ReferenceGrid->GetDims();
    const int DimsX = Dims[0], DimsY = Dims[1], DimsZ = Dims[2];

    this->Clipper.SetDeltaX( axesHashX[DimsX-1] - axesHashX[0] );
    this->Clipper.SetDeltaY( axesHashY[DimsY-1] - axesHashY[0] );
    this->Clipper.SetDeltaZ( axesHashZ[DimsZ-1] - axesHashZ[0] );
    this->Clipper.SetClippingBoundaries( this->m_FloatingCropRegionFractional );
    
    DataGrid::IndexType::ValueType startZ, endZ;
    if ( this->ClipZ( this->Clipper, axesHashZ[0], startZ, endZ ) ) 
      {
      startZ = std::max<DataGrid::IndexType::ValueType>( startZ, this->m_ReferenceCropRegion.From()[2] );
      endZ = std::min<DataGrid::IndexType::ValueType>( endZ, this->m_ReferenceCropRegion.To()[2] + 1 );
      
      const int numberOfTasks = std::min<size_t>( 4 * this->m_NumberOfThreads - 3, endZ - startZ + 1 );
      this->m_EvaluateTaskInfo.resize( numberOfTasks );
      
      for ( int taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
	{
	this->m_EvaluateTaskInfo[taskIdx].thisObject = this;
	this->m_EvaluateTaskInfo[taskIdx].AxesHash = &axesHash;
	this->m_EvaluateTaskInfo[taskIdx].StartZ = startZ;
	this->m_EvaluateTaskInfo[taskIdx].EndZ = endZ;
	}

      ThreadPool::GetGlobalThreadPool().Run( EvaluateThread, this->m_EvaluateTaskInfo );
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

  /// Metric objects for the separate threads.
  std::vector<VM> m_ThreadMetric;

  /// Mutex lock for access to global Metric field.
  MutexLock m_MetricMutex;

  /** Thread parameter block for incremental gradient computation.
   * This structure holds all thread-specific information. A pointer to an
   * instance of this structure is given to EvaluateGradientThread() for
   * each thread created.
   */
  typedef struct 
  {
    /// Pointer to the functional object that created the thread.
    Self *thisObject;
    /// Axes hash.
    const TransformedVolumeAxes* AxesHash;
    /// First plane of clipped reference volume.
    DataGrid::IndexType::ValueType StartZ;
    /// Last plane of clipped reference volume.
    DataGrid::IndexType::ValueType EndZ;
  } EvaluateTaskInfo;
 
  /// Info blocks for parallel threads evaluating functional gradient.
  std::vector<typename Self::EvaluateTaskInfo> m_EvaluateTaskInfo;

  /** Compute functional gradient as a thread.
    * This function (i.e., each thread) iterates over all parameters of the
    * current warp transformation. Among all active (i.e., not disabled)
    * parameters, it selects the ones that have an index with modulus
    * equal to the threads index when divided by the total number of threads.
    * For these parameters, the thread computes the partial derivative of the
    * functional by finite-difference approximation.
    */
  static void EvaluateThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t ) 
  {
  typename Self::EvaluateTaskInfo *info = static_cast<typename Self::EvaluateTaskInfo*>( args );

    Self *me = info->thisObject;
    const VM* Metric = me->Metric;

    VM& threadMetric = me->m_ThreadMetric[threadIdx];
    threadMetric.Reset();

    const Vector3D *hashX = (*info->AxesHash)[0], *hashY = (*info->AxesHash)[1], *hashZ = (*info->AxesHash)[2];
    Vector3D pFloating;

    const DataGrid::IndexType& Dims = me->ReferenceGrid->GetDims();
    const int DimsX = Dims[0], DimsY = Dims[1];

    int fltIdx[3];
    Types::Coordinate fltFrac[3];

    const int FltDimsX = me->FloatingDims[0], FltDimsY = me->FloatingDims[1];

    Vector3D rowStart;
    Vector3D planeStart;

    int offset;
    DataGrid::IndexType::ValueType pX, pY, pZ;
    // Loop over all remaining planes
    for ( pZ = info->StartZ + taskIdx; pZ < info->EndZ; pZ += taskCnt ) 
      {
      // Offset of current reference voxel
      int r = pZ * DimsX * DimsY;
      
      planeStart = hashZ[pZ];
      
      DataGrid::IndexType::ValueType startY, endY;
      if ( me->ClipY( me->Clipper, planeStart, startY, endY ) ) 
	{	
	startY = std::max<DataGrid::IndexType::ValueType>( startY, me->m_ReferenceCropRegion.From()[1] );
	endY = std::min<DataGrid::IndexType::ValueType>( endY, me->m_ReferenceCropRegion.To()[1] + 1 );
	r += startY * DimsX;
	
	// Loop over all remaining rows
	for ( pY = startY; pY<endY; ++pY ) 
	  {
	  (rowStart = planeStart) += hashY[pY];
	  
	  DataGrid::IndexType::ValueType startX, endX;
	  if ( me->ClipX( me->Clipper, rowStart, startX, endX ) ) 
	    {	    
	    startX = std::max<DataGrid::IndexType::ValueType>( startX, me->m_ReferenceCropRegion.From()[0] );
	    endX = std::min<DataGrid::IndexType::ValueType>( endX, me->m_ReferenceCropRegion.To()[0] + 1 );
	    
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
		threadMetric.Increment( Metric->GetSampleX( r ), Metric->GetSampleY( offset, fltFrac ) );
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

    me->m_MetricMutex.Lock();
    me->Metric->AddMetric( threadMetric );
    me->m_MetricMutex.Unlock();
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVoxelMatchingAffineFunctionalTemplate_h_included_
