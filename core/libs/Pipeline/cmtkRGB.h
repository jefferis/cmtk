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

#ifndef __cmtkRGB_h_included_
#define __cmtkRGB_h_included_

#include <cmtkTypes.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Red, green, and blue components of one pixel.
 */
typedef struct {
#ifdef WORDS_BIGENDIAN
  /// Red
  byte R;
  /// Green
  byte G;
  /// Blue
  byte B;
#else
  byte B;
  byte G;
  byte R;
#endif
} RGB;

/** RGB components plus transparency (alpha value).
 */
class RGBA {
public:
#ifdef WORDS_BIGENDIAN
  /** Transparency.
   * Opacity is linear between 0 (fully transparent) and 255 (fully opaque).
   */
  byte Alpha;
  /// Red
  byte R;
  /// Green
  byte G;
  /// Blue
  byte B;
#else
  byte B;
  byte G;
  byte R;
  byte Alpha;
#endif
  /// Assignment operator.
  void operator =( const RGB& rgb ) { R = rgb.R; G = rgb.G; B = rgb.B; }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRGB_h_included_
