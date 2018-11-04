/*
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkUniformVolumeHoughTransformSphere_h_included_
#define __cmtkUniformVolumeHoughTransformSphere_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUnits.h>

namespace cmtk {

/** \addtogroup Base */
//@{

/** Spherical Hough transform for 3D uniform image data.
 */
class UniformVolumeHoughTransformSphere {
 public:
  /// This class.
  typedef UniformVolumeHoughTransformSphere Self;

  /// Constructor: link to UniformVolume object.
  explicit UniformVolumeHoughTransformSphere(
      UniformVolume::SmartConstPtr &volume)
      : m_UniformVolume(volume) {}

  /// Get Hough-transformed data for a specific radius in world units.
  TypedArray::SmartPtr Get(
      const Types::Coordinate radius /*!< Radius of detected spheres */) const;

 private:
  /// The UniformVolume object we're working on.
  UniformVolume::SmartConstPtr m_UniformVolume;
};

}  // namespace cmtk

#endif  // #ifndef __cmtkUniformVolumeHoughTransformSphere_h_included_
