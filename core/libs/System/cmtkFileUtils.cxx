/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkFileUtils.h"

#include <limits.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#ifdef _MSC_VER
#  include <direct.h>
#  include <windows.h>
#endif

#include <iostream>

namespace 
cmtk
{

/** \addtogroup System */
//@{

namespace
FileUtils
{

int 
RecursiveMkDir( const std::string& filename, const int permissions )
{
  const int result = RecursiveMkPrefixDir( filename, permissions );
  if ( result) 
    return result;

#ifdef _MSC_VER
  return _mkdir( filename.c_str() );
#else
  return mkdir( filename.c_str(), permissions );
#endif
}

int
RecursiveMkPrefixDir
( const std::string& filename, const int permissions )
{
  char prefix[PATH_MAX];
  for ( unsigned i=0; filename[i]; ++i ) 
    {
    if ( (filename[i] == CMTK_PATH_SEPARATOR) || (filename[i] == '/') ) 
      {
      prefix[i+1] = 0;
      if ( i ) // do not delete single "/" or "\"
	prefix[i] = 0;
      else
	prefix[i] = CMTK_PATH_SEPARATOR;
      
#ifdef _MSC_VER
      if ( (i > 0) && (prefix[i-1] == ':') )
	{
	prefix[i] = '\\';
	}
#endif
#ifdef _MSC_VER
      const int result = _mkdir( prefix );
#else
      const int result = mkdir( prefix, permissions );
#endif
      if ( result && errno != EEXIST ) 
	{
	return result;
	}
      }
    prefix[i] = filename[i];
    }
  return 0;
}

std::string
GetAbsolutePath( const std::string& relPath )
{
#ifdef _MSC_VER
  char absPath[PATH_MAX];
  GetFullPathName( relPath.c_str(), PATH_MAX, absPath, NULL );
  return std::string( absPath );
#else
  if ( relPath[0] == CMTK_PATH_SEPARATOR )
    {
    return relPath;
    }
  else
    {
    char absPath[PATH_MAX];
    getcwd( absPath, PATH_MAX );
    if ( absPath[ strlen( absPath )-1 ] != CMTK_PATH_SEPARATOR )
      strcat( absPath, CMTK_PATH_SEPARATOR_STR );
    
    return std::string( absPath ) + relPath;
    }
#endif
}

} // namespace FileUtils

} // namespace cmtk
