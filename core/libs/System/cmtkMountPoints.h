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

#ifndef __cmtkMountPoints_h_included_
#define __cmtkMountPoints_h_included_

#include <cmtkconfig.h>

#include <limits.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Name of the shell variable defining the directory translation.
 * This variable can be set to contain a list of substitution rules. Every
 * ruls has the form "target=source", where all appearances of "source" are
 * to be replaced by "target" in all filesystem paths. Several of these rules
 * may be concatenated by "," -- an example is "j:=/cdrom,k:=/home". This 
 * results in all paths relative to "/cdrom" being relative to "j:" after
 * substituion. The same holds for "k:" and "/home". This example could be
 * used to read data on a PC mounting Unix filesystems to network drives.
 * Conversely, the Unix box would define "/cdrom=j:,/home=k:" in order to be
 * able to read data written by this PC.
 *@see MountPoints
 */
const char* const CMTK_MOUNTPOINTSVAR = "CMTK_MOUNTPOINTS";

/** Legacy environment variable.
 * This is what CMTK_MOUNTPOINTS used to be called in an earlier life. We
 * continue to check for it so we don't break older scripts.
 */
const char* const IGS_MOUNTPOINTSVAR = "IGS_MOUNTPOINTS";

/** Directory translation.
 * This class implements a translation for file system paths. This is used to
 * dynamically change references to the file system when using data from one
 * system on another without changing the actual files. This is useful, for
 * example, when reading cross-referenced files from a CDROM on several
 * machines with different CD mount points. It is also useful for exchanging
 * files containing paths between Unix and Windows systems.
 */
class MountPoints 
{
public:
  /** Perform directory substitutions.
   *@param The original path before substitions.
   *@return A pointer to a static buffer holding the path after all substitions
   * have been done. The buffer is guaranteed to remain unchanged until and 
   * only until the next time Translate() is called.
   *@see CMTK_MOUNTPOINTSVAR
   */
  static const char* Translate ( const char* path );

private:
  /// Static buffer holding paths after pattern substitution.
  static char Buffer[PATH_MAX];
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMountPoints_h_included_
