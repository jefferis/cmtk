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

#ifndef __cmtkViewerAnnotations_h_included_
#define __cmtkViewerAnnotations_h_included_

#include <cmtkconfig.h>

#include <cmtkTypes.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Modes for scale indicator.
 */
typedef enum {
  /// Do not draw a scale indicator.
  ANNOT_SCALE_NONE,
  /// Draw coordinate axes.
  ANNOT_SCALE_AXES,
  /// Draw a ruler.
  ANNOT_SCALE_RULER
} AnnotationScaleMode;

/** Parameter class for viewer annotations.
 */
class ViewerAnnotations 
{
public:
  /// Scale indicator mode.
  igsGetSetMacro(AnnotationScaleMode,ScaleMode);
  
  /// Scale indicator major tick distance.
  igsGetSetMacro(Types::Coordinate,ScaleMajorTickDistance);

  /// Scale indicator minor tick distance.
  igsGetSetMacro(Types::Coordinate,ScaleMinorTickDistance);

  /// Label scale indicator.
  igsGetSetMacro(bool,ScaleLabels);

  /// Flag for slice index display.
  igsGetSetMacro(bool,ShowSliceIndex);

  /// Flag for slice location display.
  igsGetSetMacro(bool,ShowSliceLocation);

  /// Constructor.
  ViewerAnnotations() 
  {
    ScaleMode = ANNOT_SCALE_NONE;
    ScaleMajorTickDistance = 5;
    ScaleMinorTickDistance = 1;
    ScaleLabels = true;
  }

  /// Copy operator.
  ViewerAnnotations& operator= ( const ViewerAnnotations& src ) 
  {
    ScaleMode = src.ScaleMode;
    ScaleMajorTickDistance = src.ScaleMajorTickDistance;
    ScaleMinorTickDistance = src.ScaleMinorTickDistance;
    ScaleLabels = src.ScaleLabels;
    return *this;
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkViewerAnnotations_h_included_
