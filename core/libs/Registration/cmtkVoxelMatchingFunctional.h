/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkVoxelMatchingFunctional_h_included_
#define __cmtkVoxelMatchingFunctional_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkFunctional.h>
#include <Base/cmtkVector.h>
#include <Base/cmtkVolume.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkMatchedLandmarkList.h>

#include <System/cmtkException.h>

#if defined(CMTK_USE_SMP)
#  include <System/cmtkLockingPtr.h>
#endif

#include <cassert>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for voxel matching functionals.
 * This class is used as the common base class for more specific functional
 * classes. It contains all data structure and functions that do not depend on
 * any template parameters introduced later in the inheritance hierarchy. It
 * should therefore help avoiding unnecessary code duplication.
 */
class VoxelMatchingFunctional :
    public Functional,
    private CannotBeCopied
{
public:
  /// This class.
  typedef VoxelMatchingFunctional Self;

  /// Superclass.
  typedef Functional Superclass;

protected:
  /// Pointer to the reference grid.
  UniformVolume::SmartPtr ReferenceGrid;

  /// Pointer to the floating grid.
  UniformVolume::SmartPtr FloatingGrid;

  /// Data class of reference image.
  DataClass ReferenceDataClass;

  /// Data class of floating image.
  DataClass FloatingDataClass;

  /// Crop region in the reference volume.
  DataGrid::RegionType m_ReferenceCropRegion;

  /// Optional list of matched landmarks.
  cmtkGetSetMacro(MatchedLandmarkList::SmartPtr,MatchedLandmarkList);

  /// Weight for the landmark registration error relative to image similarity.
  cmtkGetSetMacro(Self::ReturnType,LandmarkErrorWeight);

public:
  /** Constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   *\param reference The reference (i.e. static) volume.
   *\param floating The floating (i.e. transformed) volume.
   */
  VoxelMatchingFunctional( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
    : m_MatchedLandmarkList( NULL )
  {
    this->InitFloating( floating );
    this->InitReference( reference );
    this->m_LandmarkErrorWeight = 0;
  }

  /** Destructor.
   */
  virtual ~VoxelMatchingFunctional() {}

protected:
  /// Grid dimensions of the floating volume.
  DataGrid::IndexType FloatingDims;

  /// Extents of the floating volume in real-world coordinates.
  UniformVolume::CoordinateVectorType FloatingSize;

  /// Inverse pixel sizes of the floating volume.
  Vector3D FloatingInverseDelta;

  /// Coordinates of the floating image's cropping region.
  UniformVolume::CoordinateRegionType m_FloatingCropRegionCoordinates;

  /// Fractional index coordinates of the floating image's cropping region.
  UniformVolume::CoordinateRegionType m_FloatingCropRegionFractional;
 
  /// Grid dimensions of the reference volume.
  DataGrid::IndexType ReferenceDims;

  /// Extents of the reference volume in real-world coordinates.
  UniformVolume::CoordinateVectorType ReferenceSize;

  /// Inverse pixel deltas of the reference volume.
  UniformVolume::CoordinateVectorType ReferenceInvDelta;

  /** Find rectilinear area in original reference grid.
   *\param fromVOI Lower corner of area to find.
   *\param toVOI Upper corner of area to find.
   *\return The smallest box of reference 
   * grid voxels that contains the given coordinate range.
   */
  const DataGrid::RegionType GetReferenceGridRange ( const Vector3D& fromVOI, const Vector3D& toVOI );

private:
  /// Initialize internal data structures for floating image.
  void InitFloating( UniformVolume::SmartPtr& floating );

  /// Initialize internal data structures for reference image.
  void InitReference( UniformVolume::SmartPtr& reference );
};

/** Functional that evaluates a voxel-based similarity measure.
 * This class defines the type of functional that is optimized during
 * voxel-based registration. It holds references to reference and floating data
 * and computes similarity as well as its gradient w.r.t. a given
 * transformation.
 *
 * The metric to be optimized is given by a template parameter, therefore 
 * allowing inlined code to be generated for efficient evaluation.
 *
 * This class, however, is still universal with respect to the registration
 * transformation. Derived classes can therefore efficiently implement all
 * transformation-dependent operations.
 */
template<class VM>
class VoxelMatchingFunctional_Template 
{
protected:
  /// Pointer to the voxel-metric.
  SmartPointer<VM> Metric;

public:
  /** Constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   *\param reference The reference (i.e. static) volume.
   *\param floating The floating (i.e. transformed) volume.
   */
  VoxelMatchingFunctional_Template
  ( UniformVolume::SmartPtr& reference, 
    UniformVolume::SmartPtr& floating )
  { Metric = SmartPointer<VM>( new VM( reference, floating ) ); }
 
  /** Destructor.
   * Delete metric object.
   */
  virtual ~VoxelMatchingFunctional_Template () {}

#if defined(CMTK_USE_SMP)
protected:
   /** Mutex lock.
    * This mutex is used by the multi-threaded complete functional evaluation.
    * At the end of its computation, each thread adds its contribution in the
    * form of a local 2-D histogram to the global histogram. This one step is
    * protected by the mutex lock, as it requires exclusive write access to the
    * global metric object. It makes sense, however, to have the threads do
    * this final step as some of them may finish before others, which gives
    * them the required exclusive time.
    */
   MutexLock MetricMutex;
#endif
};

//@}

} // namespace cmtk

#endif // __cmtkVoxelMatchingFunctional_h_included_
