/*
//
//  Copyright 2012 SRI International
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

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Class for detecting landmark locations of the Magphan EMR051 structural imaging phantom.
 */
class DetectPhantomMagphanEMR051
{
public:
  /// This class.
  typedef DetectPhantomMagphanEMR051 Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const for  this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Spatial coordinate vector.
  typedef UniformVolume::SpaceVectorType SpaceVectorType;
  
  /// Constructor: compute registrations.
  DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage );
  
  /// Get landmark coordinates.
  std::vector<UniformVolume::SpaceVectorType> GetLandmarks();

  /// Get mask of detected sphere locations.
  UniformVolume::SmartConstPtr GetDetectedSpheresMask() const
  {
    return this->m_ExcludeMask;
  }

private:
  /// Image of the phantom.
  UniformVolume::SmartConstPtr m_PhantomImage;

  /** Evolving exclusion mask.
   * When a sphere is detected, its volume is marked as off-limits herein so other spheres are not incorrectly detected in the same place.
   */
  UniformVolume::SmartPtr m_ExcludeMask;

  /** Temporary inclusion mask.
   * When we detect a sphere in a specific pre-determined area, this mask contains as non-zero the candidate pixels.
   */
  UniformVolume::SmartPtr m_IncludeMask;

  /// Find a sphere of given radius.
  Self::SpaceVectorType FindSphere( const TypedArray& filterResponse );

  /// Find a sphere in a band of given radius.
  Self::SpaceVectorType FindSphereAtDistance( const TypedArray& filterResponse, const Self::SpaceVectorType& bandCenter, const Types::Coordinate bandRadius, const Types::Coordinate bandWidth );

  /// Refine sphere position based on intensity-weighted center of mass.
  Self::SpaceVectorType RefineSphereLocation( const Self::SpaceVectorType& estimate, const Types::Coordinate radius, const int margin, const int label );
};

//@}

} // namespace cmtk
