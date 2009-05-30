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

#ifndef __cmtkIntersection_h_included_
#define __cmtkIntersection_h_included_

#include <cmtkVector3D.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Utility Class for Intersection Computation.
 * This class provides functions to compute intersections of volumes with
 * lines, planes, and other volumes. Line-to-volume intersection can be
 * computed by a single function call. 
 *
 * Plane to volume intersection is done by first computing a range of parallel
 * lines that intersect the volume. Then, for each line the line to volume
 * intersection has to be computed to determine the exact intersection.
 *
 * Analog, for volume to volume intersection, first a range of parallel planes
 * is determined. Those then have to be examined in separate as described 
 * before.
 *
 * This class' main application is ray clipping for DRR computation and volume
 * clipping for voxel based volume similarity computation. All member functions
 * are static, so they can be used without contructing an object.
 *@author T. Rohlfing
 *@version $Revision$ $Date$
 */
class Intersection 
{
public:
  /** Compute line-to-volume intersection.
   * This function computes the part of a lines that lies within a given
   * volume.
   *@return This function returns 1 if and only if there is a non-empty
   * intersection of line and volume. In this case, the intersection is
   * described by fromFactor and toFactor.
   *@param fromFactor If the function returned 1, this variable holds the
   * relative distance to the entrance point of the line into the volume.
   * 0 means the line's starting point, 1 means its end.
   *@param toFactor If the function returned 1, this variable holds the 
   * relative distance to the exit point of the line from the volume. Possible
   * values range from 0 to 1 and have the same meaning as fromFactor.
   *@param offset This 3D vector is the line's starting point.
   *@param dX This is the vector spanning from the starting point of the line 
   * to its end.
   *@param Size This is the size of the volume. It is always rectilinear with
   * faces parallel to the coordinate axes and its lower left front corner 
   * identical to the coordinate origin. If affine transformed volumes are
   * to be tested, the inverse affine transformation has to be applied to
   * the line's endpoints instead.
   *@param initFromFactor The fromFactor parameter's value is initialized
   * with this value. It is therefore the lower bound of the parameter range
   * that is available for intersection.
   *@param initToFactor The toFactor parameter's value is initialized
   * with this value. It is therefore the upper bound of the parameter range
   * that is available for intersection. One application for this parameter
   * is to use a value bigger than 1, even if [0,1] is the allowed range. Then,
   * by testing if toFactor == 1.0, it can be determined whether clipping
   * set the value to 1. This, for example allows to tell closed from open
   * intervals, which may be important for subsequent computation such as
   * volume probing.
   *@param lowerClosed This flag defines whether lower range boundaries are
   * open (value 0, default) or closed (value 1). In case of an open range, 
   * the bounding value itself is not an element of the range. Thus, if an
   * intersection is entirely on the boundary, then it is empty in case of an
   * open range.
   *@param upperClosed This flag has the same meaning and default as
   * lowerClosed, but refers to the ranges' upper bounds.
   */
  static int IntersectX ( Types::Coordinate& fromFactor, Types::Coordinate& toFactor,
			  const Vector3D& offset, const Vector3D& dX, 
			  const Types::Coordinate Size[3],
			  const Types::Coordinate initFromFactor = 0,
			  const Types::Coordinate initToFactor = 1,
			  const int lowerClosed = 0,
			  const int upperClosed = 0 );

  /** Compute plane-to-volume intersection.
   * This function computes the part of a plane that intersects with a given
   * volume. The intersection is only tested with respect to the dY direction.
   * for the range of relative positions returned by fromFactor and toFactor,
   * IntersectX can be used to computed the exact 2D intersection.
   *
   * Parameters and return value are identical to IntersectionX.
   *@param dY This is the second vector in addition to dX spanning the plane
   * under consideration. In general, dX and dY should be orthogonal, and dY
   * should be defined in the way just as dX is, i.e. as the difference of two
   * plane corner vectors.
   */
  static int IntersectY ( Types::Coordinate& fromFactor, Types::Coordinate& toFactor,
			  const Vector3D& offset, const Vector3D& dX, 
			  const Vector3D& dY, const Types::Coordinate Size[3],
			  const Types::Coordinate initFromFactor = 0,
			  const Types::Coordinate initToFactor = 1 );

  /** Compute volume-to-volume intersection.
   * This function computes the part of a volume that intersects with a given
   * volume. The intersection is only tested with respect to the dZ direction.
   * for the range of relative positions returned by fromFactor and toFactor,
   * IntersectY and afterwards IntersectX can be used to computed the exact 3D
   * intersection.
   *
   * Parameters and return value are identical to IntersectionX.
   *@param dY This is the third vector in addition to dX and dY spanning the
   * volume under consideration. In general, dX, dY, and dZ should be
   * orthogonal, and dZ should be defined in the same way as dX and dY are, 
   * i.e. as the difference of two plane corner vectors where one is the volume
   * origin.
   */
  static int IntersectZ ( Types::Coordinate& fromFactor, Types::Coordinate& toFactor,
			  const Vector3D& offset, const Vector3D& dX, 
			  const Vector3D& dY, const Vector3D& dZ,
			  const Types::Coordinate Size[3],
			  const Types::Coordinate initFromFactor = 0,
			  const Types::Coordinate initToFactor = 1 );

};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkIntersection_h_included_
