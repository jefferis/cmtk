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

#ifndef __cmtkPlane_h_included_
#define __cmtkPlane_h_included_

#include <cmtkconfig.h>

#include <Pipeline/cmtkPipelineObject.h>

#include <Base/cmtkMacros.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkVector3D.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class for 2D planes, that is uniform point meshes.
 */
class Plane : 
  public PipelineObject
{
public:
  /// Create new object.
  static Plane* New() { return new Plane; }

  /// Dimensions array.
  igsClassParameter2Array(unsigned int,Dims);

  /// Spacing (ie. pixel size) array.
  igsClassParameter2Array(Types::Coordinate,Spacing);

  /// Origin of image in 3D space.
  igsClassParameter3Array(Types::Coordinate,Origin);

  /// Direction of image's x-axis in 3D space.
  igsClassParameter3Array(Types::Coordinate,DirectionX);

  /// Direction of image's y-axis in 3D space.
  igsClassParameter3Array(Types::Coordinate,DirectionY);

  /** Copy the structure of another Plane object.
   * This function copies dimensions, pixel size, and spatial location of
   * a given object.
   */
  void CopyStructure( const Plane *plane );

  /// Return number of pixels in this object.
  virtual unsigned int GetNumPixels() const { return Dims[0]*Dims[1]; }

  /// Return 3D coordinate of a particular pixel.
  void GetPixelLocation( Vector3D& v, const unsigned int x, const unsigned int y ) const 
  {
    for ( int dim = 0; dim<3; ++dim )
      v[dim] = Origin[dim] + x * DirectionX[dim] * Spacing[0] + y * DirectionY[dim] * Spacing[1];
  }

  /** Project 3D coordinate onto plane.
   *@param p Projected coordinate.
   *@param q Original coordinate.
   */
  void Project( Vector3D& p, const Vector3D& q ) const;
  
  /** Project 3D coordinate onto image plane pixels.
   *@param v Original coordinate.
   *@param i Index of projected pixel in x direction.
   *@param j Index of projected pixel in y direction.
   */
  void ProjectPixel( const Vector3D& v, unsigned int& i, unsigned int& j ) const;
  
protected:
  /** Default constructor.
   * Set all plane fields to safe values.
   */
  Plane();

  /// Destructor.
  virtual ~Plane() {};
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkPlane_h_included_
