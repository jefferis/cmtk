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

#include "cmtkMountPoints.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

char MountPoints::Buffer[PATH_MAX];

const char* 
MountPoints::Translate( const char* path )
{
  // Get environment variable defining mount points.
  const char *mountpoints = getenv( CMTK_MOUNTPOINTSVAR );
  if ( ! mountpoints )
    mountpoints = getenv( IGS_MOUNTPOINTSVAR );

  // Not defined: Return path unmodified
  if ( ! mountpoints ) return path;
  strcpy( Buffer, path );

  char searchStr[256], replaceStr[256];
  const char *delim;

  const char *nextRule = mountpoints;

  while ( nextRule ) 
    {
    delim = strchr( nextRule, '=' );
    if ( delim ) 
      {
      int cplen = delim - nextRule;
      strncpy( searchStr, nextRule, cplen );
      searchStr[cplen] = 0;
      
      nextRule = strchr( delim, ',' );
      if ( nextRule ) 
	{
	int cplen = nextRule - delim - 1;
	strncpy( replaceStr, delim+1, cplen );
	replaceStr[cplen] = 0;
	nextRule++;
	} 
      else
	{
	strcpy( replaceStr, delim+1 );
	nextRule = NULL;
	}
      
      // check for beginning-of-line token
      bool checkPrefixOnly = false;
      if ( searchStr[0] == '^' ) 
	{
	checkPrefixOnly = true;
	}
      
      if ( checkPrefixOnly ) 
	{
	// Check if rule applies to given path.
	if ( !strncmp( path, searchStr+1, strlen( searchStr ) - 1 ) ) 
	  {
	  // Yes, it does: Substitute prefix accordingly and return pointer
	  // to buffer containing modified path.
	  strcat( strcpy( Buffer, replaceStr ), path+strlen(searchStr)-1 );
	  return Buffer;
	  }
	} 
      else
	{
	// Substitute non-prefix occurences as well
	char *found = NULL;
	if ( ( found = strstr( Buffer, searchStr ) ) ) 
	  {
	  // Yes, it does: Substitute accordingly and return pointer
	  // to buffer containing modified path.
	  char tmpPath[PATH_MAX];
	  memset( tmpPath, 0, sizeof( tmpPath ) );
	  strcat( strcat( strncpy( tmpPath, Buffer, found-Buffer ), replaceStr ), found + strlen(searchStr) );
	  strcpy( Buffer, tmpPath );
	  }
	}
      }
    }
  
  return Buffer;
}

} // namespace cmtk
