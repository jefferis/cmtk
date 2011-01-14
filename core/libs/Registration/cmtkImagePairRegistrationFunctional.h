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

#ifndef __cmtkImagePairRegistrationFunctional_h_included_
#define __cmtkImagePairRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkFunctional.h>
#include <Base/cmtkVector.h>
#include <Base/cmtkVolume.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkMatchedLandmarkList.h>

#include <System/cmtkException.h>
#include <System/cmtkLockingPtr.h>

#include <Registration/cmtkImagePairSimilarityMeasure.h>

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
class ImagePairRegistrationFunctional :
  /// Inherit Functional interface.
  public Functional 
{
public:
  /// This class.
  typedef ImagePairRegistrationFunctional Self;

  /// Superclass.
  typedef Functional Superclass;

protected:
  /// Pointer to the reference grid.
  UniformVolume::SmartConstPtr m_ReferenceGrid;

  /// Pointer to the floating grid.
  UniformVolume::SmartConstPtr m_FloatingGrid;

  /// Data class of reference image.
  DataClass m_ReferenceDataClass;

  /// Data class of floating image.
  DataClass m_FloatingDataClass;

  /// Rectangular crop region in the reference volume.
  DataGrid::RegionType m_ReferenceCropRegion;

  /// Optional list of matched landmarks.
  cmtkGetSetMacro(MatchedLandmarkList::SmartPtr,MatchedLandmarkList);

  /// Weight for the landmark registration error relative to image similarity.
  cmtkGetSetMacro(Self::ReturnType,LandmarkErrorWeight);

public:
  /** Constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   *@param reference The reference (i.e. static) volume.
   *@param floating The floating (i.e. transformed) volume.
   */
  ImagePairRegistrationFunctional( UniformVolume::SmartConstPtr& reference, UniformVolume::SmartConstPtr& floating )
    : m_MatchedLandmarkList( NULL ),
      m_ForceOutsideFlag( false )
  {
    this->InitFloating( floating );
    this->InitReference( reference );
    this->m_LandmarkErrorWeight = 0;
  }

  /** Destructor.
   */
  virtual ~ImagePairRegistrationFunctional() {}

  /// Set flag and value for forcing values outside the floating image.
  /// Set flag and value for forcing values outside the floating image.
  virtual void SetForceOutside
  ( const bool flag = true, const Types::DataItem value = 0 )
  {
    this->m_ForceOutsideFlag = flag;
    this->m_ForceOutsideValueRescaled = this->m_Metric->GetFloatingValueScaled( value );
  }

protected:
  /// The metric (similarity measure) object.
  ImagePairSimilarityMeasure::SmartPtr m_Metric;

  /// Grid dimensions of the floating volume.
  DataGrid::IndexType m_FloatingDims;

  /// Extents of the floating volume in real-world coordinates.
  UniformVolume::CoordinateVectorType m_FloatingSize;

  /// Inverse pixel sizes of the floating volume.
  UniformVolume::CoordinateVectorType m_FloatingInverseDelta;

  /// Coordinates of the floating image's cropping region.
  UniformVolume::CoordinateRegionType m_FloatingCropRegionCoordinates;
 
  /// Fractional index starting coordinate of the floating's cropping region.
  UniformVolume::CoordinateRegionType m_FloatingCropRegionFractIndex;

  /// Grid dimensions of the reference volume.
  DataGrid::IndexType m_ReferenceDims;

  /// Extents of the reference volume in real-world coordinates.
  UniformVolume::CoordinateVectorType m_ReferenceSize;

  /// Inverse pixel deltas of the reference volume.
  UniformVolume::CoordinateVectorType m_ReferenceInverseDelta;

  /// Flag for forcing pixel values outside the floating image.
  bool m_ForceOutsideFlag;

  /// Rescaled byte value for forcing pixel values outside the floating image.
  Types::DataItem m_ForceOutsideValueRescaled;

  /** Find rectilinear area in original reference grid.
   *@param fromVOI Lower corner of area to find.
   *@param toVOI Upper corner of area to find.
   *@return The smallest region of reference grid voxels that contains the given coordinate range.
   */
  const DataGrid::RegionType GetReferenceGridRange ( const Vector3D& fromVOI, const Vector3D& toVOI );

private:
  /// Private copy constructor: prevent copying.
  ImagePairRegistrationFunctional ( const ImagePairRegistrationFunctional& ) {}

  /// Initialize internal data structures for floating image.
  void InitFloating( UniformVolume::SmartConstPtr& floating );

  /// Initialize internal data structures for reference image.
  void InitReference( UniformVolume::SmartConstPtr& reference );
};

//@}

} // namespace cmtk

#endif // __cmtkImagePairRegistrationFunctional_h_included_
