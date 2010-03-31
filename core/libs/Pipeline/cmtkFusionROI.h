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

#ifndef __cmtkFusionROI_h_included_
#define __cmtkFusionROI_h_included_

#include "cmtkArrayFilter.h"

#include "cmtkImageRGB.h"
#include "cmtkMacros.h"
#include "cmtkRGB.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/// Class for Region-of-Interest image fusion.
class FusionROI : 
  /// This is a filter from two RGB images to an RGB image.
  public ArrayFilter<ImageRGB,ImageRGB,2> 
{
public:
  /// Create new class instance.
  static FusionROI* New() { return new FusionROI; }

  /** Parameter selecting which ROI fusion technique to use.
   * Supported values are: 0 - Vertical bars, 1 - Horizontal bars, 2 - Circle,
   * 3 - Checkerboard.
   */
  igsClassParameter(int,Mode);

  /// Parameter selecting which of both images to display "on top".
  igsClassParameter(int,TopImageIndex);

  /**@name "Vertical Bars"-mode parameters.
   */
  //@{
  /// Relative bar width.
  igsClassParameter(double,VerticalWidth);

  /// Relative bar offset
  igsClassParameter(double,VerticalOffset);
  //@}

  /**@name "Horizontal Bars"-mode parameters.
   */
  //@{
  /// Relative bar height.
  igsClassParameter(double,HorizontalHeight);

  /// Relative bar offset.
  igsClassParameter(double,HorizontalOffset);
  //@}

  /**@name "Circle"-mode parameters.
   */
  //@{
  /// Relative x-coordinate of the circle center.
  igsClassParameter(double,CircleX);

  /// Relative y-coordinate of the circle center.
  igsClassParameter(double,CircleY);

  /// Relative radius of the circle.
  igsClassParameter(double,CircleDelta);

  /// Relative offset of the circle.
  igsClassParameter(double,CircleOffset);
  //@}

  /**@name "Box"-mode parameters.
   */
  //@{
  /// Relative x-coordinate of the box center.
  igsClassParameter(double,BoxX);

  /// Relative y-coordinate of the box center.
  igsClassParameter(double,BoxY);

  /// Relative x-size of the box.
  igsClassParameter(double,BoxDeltaX);

  /// Relative y-size of the box.
  igsClassParameter(double,BoxDeltaY);
  //@}

  /**@name "Checkerboard"-mode parameters.
   */
  //@{
  /// Relative width of the checkerboard fields.
  igsClassParameter(double,CheckerWidth);

  /// Relative height of the checkerboard fields.
  igsClassParameter(double,CheckerHeight);

  /// Relative offset of the checkerboard fields in x-direction.
  igsClassParameter(double,CheckerOffsetX);

  /// Relative offset of the checkerboard fields in y-direction.
  igsClassParameter(double,CheckerOffsetY);
  //@}

  /// Perform fusion.
  virtual void Execute();

protected:
  /// Constructor.
  FusionROI();

  /// Virtual destructor.
  virtual ~FusionROI() {};

private:
  /// Auxiliary function to generate multiple ROI circles.
  void Circles ( RGB *target, const RGB *source1, const RGB *source2, const int nx, const int ny );

  /// Auxiliary function to generate a single ROI circle.
  void Circle ( RGB *target, const RGB *source1, const RGB *source2, const int nx, const int ny );

  /// Auxiliary function to generate a single ROI box.
  void Box ( RGB *target, const RGB *source1, const RGB *source2, const int nx, const int ny );

  /// Auxiliary function to generate a vertical bar ROI pattern.
  void VerticalSlice ( RGB* target, const RGB *source1, const RGB *source2, const int nx, const int ny, const int xwidth = -1, const int xoffset = -1 );

  /// Auxiliary function to generate a horizontal bar ROI pattern.
  void HorizontalSlice ( RGB* target, const RGB *source1, const RGB *source2, const int nx, const int ny );

  /// Auxiliary function to generate a checkerboard ROI pattern.
  void Checkerboard ( RGB* target, const RGB *source1, const RGB *source2, const int nx, const int ny );
};

#endif // #ifndef __cmtkFusionROI_h_included_

} // namespace cmtk
