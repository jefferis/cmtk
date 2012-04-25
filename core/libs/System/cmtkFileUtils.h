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

#ifndef __cmtkFileUtil_h_included_
#define __cmtkFileUtil_h_included_

#include <cmtkconfig.h>

#include <string>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Utility functions for file and directory access.
 */
namespace FileUtils
{
/** Recursively create a new directory.
 *\return Zero if operation was successful, mkdir() error code otherwise.
 */
int RecursiveMkDir( const std::string& filename, const int permissions = 0755 );

/** Recursively create all directories on a full filename's prefix.
 *\return Zero if operation was successful, mkdir() error code otherwise.
 */
int RecursiveMkPrefixDir( const std::string& filename, const int permissions = 0755 );

/// Make an absolute path name from a (possibly) relative path.
std::string GetAbsolutePath( const std::string& relPath );
}

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFileUtil_h_included_
