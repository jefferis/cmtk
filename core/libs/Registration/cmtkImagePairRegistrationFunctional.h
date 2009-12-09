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

#ifndef __cmtkImagePairRegistrationFunctional_h_included_
#define __cmtkImagePairRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkMathUtil.h>
#include <cmtkTypes.h>

#include <cmtkFunctional.h>

#include <cmtkVector.h>
#include <cmtkRect3D.h>
#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>
#include <cmtkMatchedLandmarkList.h>

#include <cmtkException.h>

#if defined(CMTK_BUILD_SMP)
#  include <cmtkLockingPtr.h>
#endif

#include <assert.h>

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
  UniformVolume::SmartPtr ReferenceGrid;

  /// Pointer to the floating grid.
  UniformVolume::SmartPtr FloatingGrid;

  /// Data class of reference image.
  DataClass ReferenceDataClass;

  /// Data class of floating image.
  DataClass FloatingDataClass;

  /// Beginning of the rectangular crop region in the reference volume.
  int ReferenceCropFrom[3];

  /// End of the rectangular crop region in the reference volume.
  int ReferenceCropTo[3];

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
  ImagePairRegistrationFunctional( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
    : m_MatchedLandmarkList( NULL )
  {
    this->InitFloating( floating );
    this->InitReference( reference );
    this->m_LandmarkErrorWeight = 0;
  }

  /** Copy constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   */
  ImagePairRegistrationFunctional ( ImagePairRegistrationFunctional& source ) : 
    m_MatchedLandmarkList( source.m_MatchedLandmarkList )
  {
    this->InitFloating( source.FloatingGrid );
    this->InitReference( source.ReferenceGrid );
    this->m_LandmarkErrorWeight = source.m_LandmarkErrorWeight;
  }

  /** Copy constructor.
   * Init pointers to volume and transformation objects and initialize
   * internal data structures.
   */
  ImagePairRegistrationFunctional ( ImagePairRegistrationFunctional *const source ) : 
    m_MatchedLandmarkList( source->m_MatchedLandmarkList )
  {
    this->InitFloating( source->FloatingGrid );
    this->InitReference( source->ReferenceGrid );
    this->m_LandmarkErrorWeight = source->m_LandmarkErrorWeight;
  }

  /** Destructor.
   */
  virtual ~ImagePairRegistrationFunctional() {}

protected:
  /// Grid dimensions of the floating volume.
  int FloatingDims[3];

  /// Extents of the floating volume in real-world coordinates.
  Types::Coordinate FloatingSize[3];

  /// Inverse pixel sizes of the floating volume.
  Vector3D FloatingInverseDelta;

  /// Starting coordinate of the floating's cropping region.
  Types::Coordinate FloatingCropFrom[3];

  /// End coordinate of the floating's cropping region.
  Types::Coordinate FloatingCropTo[3];
 
  /// Fractional index starting coordinate of the floating's cropping region.
  Types::Coordinate FloatingCropFromIndex[3];

  /// Fractional index end coordinate of the floating's cropping region.
  Types::Coordinate FloatingCropToIndex[3];
 
  /// Grid dimensions of the reference volume.
  int ReferenceDims[3];

  /// Extents of the reference volume in real-world coordinates.
  Types::Coordinate ReferenceSize[3];

  /// Inverse pixel deltas of the reference volume.
  Types::Coordinate ReferenceInvDelta[3];

  /** Find rectilinear area in original reference grid.
   *@param fromVOI Lower corner of area to find.
   *@param toVOI Upper corner of area to find.
   *@param startX On return, this reference holds the index in x direction of
   * the original reference grid that is the LOWER bound of the region defined
   * by fromVOI and toVOI. The parameters startY and startZ have equivalent
   * meanings.
   *@param voi On return, this reference holds the smallest box of reference 
   * grid voxels that contains the given coordinate range.
   */
  void GetReferenceGridRange ( const Vector3D& fromVOI, 
			       const Vector3D& toVOI,
			       Rect3D& voi );

private:
  /// Initialize internal data structures for floating image.
  void InitFloating( UniformVolume::SmartPtr& floating );

  /// Initialize internal data structures for reference image.
  void InitReference( UniformVolume::SmartPtr& reference );
};

//@}

} // namespace cmtk

#endif // __cmtkImagePairRegistrationFunctional_h_included_
