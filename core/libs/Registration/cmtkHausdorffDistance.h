/*
//
//  Copyright 2011 SRI International
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

#ifndef __cmtkHausdorffDistance_h_included_
#define __cmtkHausdorffDistance_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for computing Hausdorff distance between two label images.
 * Distance computation is implemented via the Euclidean distance maps of the two images.
 */
class HausdorffDistance
{
public:
  /// This class.
  typedef HausdorffDistance Self;

  /// Constructor.
  HausdorffDistance( UniformVolume::SmartConstPtr& image0, UniformVolume::SmartConstPtr& image1 );

  /// Get distance of two binary label maps.
  Types::Coordinate GetBinary() const;

private:
  /// First image.
  UniformVolume::SmartConstPtr m_Image0;

  /// Second image.
  UniformVolume::SmartConstPtr m_Image1;

  /// Utility function: compute "half" (i.e., one direction) of the distance term from an image (treated as binary map) and the distance map of the other image.
  static Types::Coordinate HalfDistanceBinary( const UniformVolume& image, const UniformVolume& dmap );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkHausdorffDistance_h_included_
