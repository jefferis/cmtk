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

#ifndef __cmtkPlaneSource_h_included_
#define __cmtkPlaneSource_h_included_

#include <cmtkconfig.h>

#include <cmtkMultiFilter.h>

#include <cmtkPlane.h>
#include <cmtkVolumeWrapper.h>

#include <cmtkTypes.h>
#include <cmtkMacros.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class to generate slicing planes according to volume data.
 */
class PlaneSource : 
  public MultiFilter<Plane> 
{
public:
//@{

  typedef enum
  {
/// Axial scans from foot to head.
    SCANDIRECTION_CAUDAL_CRANIAL = 0,
/// Axial scans from head to foot.
    SCANDIRECTION_CRANIAL_CAUDAL = 1,
/// Sagittal scans from right to left.
    SCANDIRECTION_RIGHT_LEFT = 2,
/// Sagittal scans from left to right.
    SCANDIRECTION_LEFT_RIGHT = 3,
/// Coronal scans from front to back.
    SCANDIRECTION_VENTRAL_DORSAL = 4,
/// Coronal scans from back to front.
    SCANDIRECTION_DORSAL_VENTRAL = 5
  } ScandirectionEnum;

  /// Create new object.
  static PlaneSource* New() { return new PlaneSource; }

  /// Return virtual class name.
  virtual const char* GetClassName() const { return "PlaneSource"; }

  /** Direction of the orthogonal reformatted plane.
   * This flag selects one of the six standard orthogonal slice plane 
   * directions. There are three pairs of directions, where each consists of
   * two identical position with a flip in the x-direction. For convenient
   * selection of the scan direction, SCANDIRECTION_XXX preprocessor definitions 
   * are available.
   *@see SCANDIRECTION_XXX
   */
  igsClassParameter(int,Direction);

  /** Position of the slice plane with respect to the selected scan direction.
   * This value is the distance of the slice plane from the reference volume 
   * origin in direction of the plane normal. The unit is identical to the 
   * volume's coordinate units, ie. it is usually in [mm].
   */
  igsClassParameter(Types::Coordinate,Position);

  /** Resolution of the generated plane.
   * This is the size of the square grid elements (pixels) of the plane object
   * created by this class.
   */
  igsClassParameter(Types::Coordinate,Resolution);

  /** Flag selection reference volume.
   * The resulting plane will be aligned with the boundaries of the given 
   * volume in the Input[] list.
   */
  igsClassParameter(int,ReferenceVolumeIndex);

  /** Execute plane generation.
   */
  virtual void Execute();

  /// Return the minimum allowed position value.
  Types::Coordinate GetMinPosition();

  /// Return the maximum allowed position value.
  Types::Coordinate GetMaxPosition();

  /// Return maximum resolution (minimum pixel size) in [mm].
  Types::Coordinate GetMaxResolution();

  /// Set input volume(s).
  virtual void SetInput( const int index, VolumeWrapper *const input );

  /// Check whether we have all valid inputs.
  virtual int HasValidInputs() const;

protected:
  /// Default constructor.
  PlaneSource();

  /// Destructor.
  virtual ~PlaneSource() {};

  /// Array of input volumes to determine the parameter ranges.
  VolumeWrapper *Input[2];
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkPlaneSource_h_included_
