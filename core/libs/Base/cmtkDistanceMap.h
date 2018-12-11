/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkDistanceMap_h_included_
#define __cmtkDistanceMap_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Base class for distance map template class.
 * This class provides typedefs and static constants that do not depend on the
 * template parameter of the derived class template.
 */
class DistanceMap
{
public:
    /// Constant used to mark unprocessed pixels.
  static const long int EDT_MAX_DISTANCE_SQUARED = 2147329548;

  /** Enumeration with binary flags that control distance map computation.
   * The defined values can be combined by arithmetic "or".
   */
  typedef enum 
  {
    /** No special functions. 
     * This flag will create a distance-from-feature map where any non-zero 
     * voxel in the feature image will be considered.
     */
    DEFAULT = 0,
    /** Compute distance-from-background map.
     * The resulting distance values will be non-zero at all feature voxels,
     * specifying the distance to the nearest non-feature voxel.
     */
    INSIDE = 1,
    /** Compute "inside", rather than "outside" map. Ths essentially inverts
     * the mask.
     */
    VALUE_EXACT = 2,
    /** Use specific feature value.
     * If this flag is set, only voxels in the feature image with values equal
     * to or larger than a given constant will be considered feature voxels.
     * All voxels with lower values will be considered background voxels.
     */
    VALUE_WINDOW = 4,
    /** Use window around specific feature value.
     * If this flag is set, only voxels in the feature image with values that are
     * within a range from a given constant will be considered feature voxels.
     * All voxels with lower values will be considered background voxels.
     */
    VALUE_THRESHOLD = 8,
    /** Compute signed distance map.
     * The "INSIDE" flag determines whether negative distance values are assigned to
     * pixels inside (flag off) or outside (flag on) the labelled region.
     */
    SIGNED = 16,
    /** Compute squared distance - do not apply final sqrt() operator.
     * This can increase efficiency if the outside code wants the squared distance in the first place.
     */
    SQUARED = 32
  } Flags;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDistanceMap_h_included_
