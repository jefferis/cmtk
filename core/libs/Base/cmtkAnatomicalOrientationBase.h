/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkAnatomicalOrientationBase_h_included_
#define __cmtkAnatomicalOrientationBase_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Base class for handling anatomical image orientation.
class AnatomicalOrientationBase
{
public:
  /// This class.
  typedef AnatomicalOrientationBase Self;

  /// Orientation of standardized reoriented image axes (LR/AP/IS).
  static const char *const ORIENTATION_STANDARD;

  /// Standard CMTK coordinate space (LR/PA/IS).
  static const char *const SPACE_CMTK;

  /// Standard ITK coordinate space (RL/AP/IS).
  static const char *const SPACE_ITK;

  /** Get closest orientation from a list.
   * This function is used to determine which orientation to bring an image into so it can be written to a file
   * format with limited orientation support (e.g., Analyze).
   */
  static const char* GetClosestOrientation( const char* desiredOrientation, const char *const availableOrientations[] );

  /** Return true if the direction corresponding to the 
   * character 'from' is on the same axis as that corresponding
   * to 'to'.
   *\param from Either L, R, A, P, I, or S
   *\param to Either L, R, A, P, I, or S 
   */
  static bool OnSameAxis( const char from, const char to );

protected:
  /// Get inverse of axis orientation.
  static char OppositeDirection( const char direction )
  {
    const char table[27] = "PbcdefghSjkRmnoAqLItuvwxyz";
    return table[direction-'A'];    
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkAnatomicalOrientationBase_h_included_
