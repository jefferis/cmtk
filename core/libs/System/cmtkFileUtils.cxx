/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <string.h>
#include <limits.h>
#include <stdlib.h>

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
#endif

namespace 
cmtk
{

/** \addtogroup System */
//@{

namespace
FileUtils
{

int 
RecursiveMkDir( const char *filename, const int permissions )
{
  int result = RecursiveMkPrefixDir( filename, permissions );
  if ( result) return result;
#ifdef _MSC_VER
  return _mkdir( filename );
#else
  return mkdir( filename, permissions );
#endif
}

int
RecursiveMkPrefixDir
( const char *filename, const int permissions )
{
  char prefix[PATH_MAX];
  struct stat buf;
  for ( unsigned i=0; filename[i]; ++i ) 
    {
    if ( (filename[i] == CMTK_PATH_SEPARATOR) ) 
      {
      prefix[i] = 0;
      if ( stat( prefix, &buf ) != 0 ) 
	{
#ifdef _MSC_VER
	if ( (i > 0) && (prefix[i-1] == ':') )
	  {
	  prefix[i] = '\\';
	  prefix[i+1] = 0;
	  }
	int result = _mkdir( prefix );
#else
	int result = mkdir( prefix, permissions );
#endif
	if ( result ) return result;
	}
      }
    prefix[i] = filename[i];
    }
  return 0;
}

char* 
GetAbsolutePath( char *absPath, const char* relPath )
{
  if ( relPath[0] == CMTK_PATH_SEPARATOR )
    {
    strcpy( absPath, relPath );
    }
  else
    {
    getcwd( absPath, PATH_MAX );
    if ( absPath[ strlen( absPath )-1 ] != CMTK_PATH_SEPARATOR )
      strcat( absPath, CMTK_PATH_SEPARATOR_STR );
    
    strcat( absPath, relPath );
    }
  
  return absPath;
}

} // namespace FileUtils

} // namespace cmtk
