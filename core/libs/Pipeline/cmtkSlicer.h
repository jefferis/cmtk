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

#ifndef __cmtkSlicer_h_included_
#define __cmtkSlicer_h_included_

#include <cmtkconfig.h>

#include <Pipeline/cmtkFilter.h>

#include <Base/cmtkMacros.h>

#include <Pipeline/cmtkPlane.h>
#include <Pipeline/cmtkImage.h>
#include <Pipeline/cmtkVolumeWrapper.h>
#include <Base/cmtkInterpolator.h>

#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Reformat slice image from a volume.
 */
class Slicer :
  /// This filter takes VolumeWrapper input and generates Image output.
  public Filter<VolumeWrapper,Image> 
{
public:
  /// Create new object.
  static Slicer* New() { return new Slicer; }

  /// Set slice plane defining geometry of reformatted data.
  void SetPlane( Plane *const plane );

  /// Flag whether or not to apply existing warp transformations.
  igsClassParameter(bool,ApplyWarp);

  /// Interpolation mode.
  igsClassParameter(cmtk::Interpolators::InterpolationEnum,InterpolationMode);

  /// Zoom factor.
  igsClassParameter(double,ZoomFactor);

  /** Check for update.
   * This functions first checks this class' additional "Plane" input object
   * for updates, then it calls the inherited "Update" function.
   *\see Object#Update
   */
  virtual long Update();

  /** Perform actual reformatting.
   */
  virtual void Execute();

protected:
  /// Default constructor.
  Slicer();

  /// Virtual destructor.
  virtual ~Slicer();

private:
  /// The plane object defining the geometry for the reformatted data.
  Plane *m_Plane;

  /// If zoom is active, this objects holds the original (unzoomed) data.
  Image *TempImage;

  /** Helper function for slicing under a linear warp transformation.
   * This function is called by Execute() if a B-spline deformation is to be
   * applied during slicing. Its main purpose is to keep separate this
   * reslicing from the purely affine implementation in Execute() itself.
   */
  void ExecuteSplineWarp( TypedArray::SmartPtr& data, const SplineWarpXform* warpXform, const unsigned int* dims, const Vector3D& offset, const Vector3D& dX, const Vector3D& dY );

  /// Convenience typedef for this class' parent.
  typedef Filter<VolumeWrapper,Image> Superclass;
};

//@}

} // namespace cmtk

#endif // __cmtkSlicer_h_included_
