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
      int result = 0;
      // NOTE(fschulze@vrvis.at): prevent to call _mkdir on drive letters 
      // (e.g. "C:\") because this would fail and subsequentually make the 
      // whole function fail Furthermore drives are not folders and cannot 
      // be created.
      
      // The current prefix describes a drive if the second letter is a colon 
      // and we only have 3 letters 
      const bool isDrive = (prefix[i-1]==':') && (i==2); 
      if (!isDrive) 
	{
        result = _mkdir( prefix );
	}
#else
      const int result = mkdir( prefix, permissions );
#endif
      if ( result && errno != EEXIST && errno != EISDIR ) 
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

std::string
Basename( const std::string& path, const std::string& suffix )
{
  std::string basename = path;

  // if given, and if present, remove suffix
  if ( ! suffix.empty() && (basename.length() >= suffix.length() ) )
    {
    if ( basename.compare( basename.length() - suffix.length(), suffix.length(), suffix ) )
      {
      basename = basename.substr( 0, basename.length() - suffix.length() );
      }
    }

  // split into dirname and basename now
  const size_t separator = basename.rfind( CMTK_PATH_SEPARATOR );
  if ( separator == std::string::npos )
    {
    // no separator, return the whole thing
    return basename;
    }
  else
    {
    // seprator - return whatever comes after
    return basename.substr( separator+1 );
    }
}

} // namespace FileUtils

} // namespace cmtk
