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

#include "cmtkMountPoints.h"

#include <stdlib.h>
#include <string.h>
#include <limits.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

std::string
MountPoints::Translate( const std::string& path )
{
  // Get environment variable defining mount points.
  const char *mountpoints = getenv( CMTK_MOUNTPOINTSVAR );
  if ( ! mountpoints )
    {
    mountpoints = getenv( IGS_MOUNTPOINTSVAR );
    
    // Not defined: Return path unmodified
    if ( ! mountpoints ) 
      return path;
    }

  std::string buffer = path;

  const char *nextRule = mountpoints;
  while ( nextRule ) 
    {
    // find the equation sign between pattern and replacement
    const char* delim = strchr( nextRule, '=' );
    if ( delim ) 
      {
      const int cplen = delim - nextRule;
      const std::string pattern = std::string( nextRule ).substr( 0, cplen );
      std::string replacement = std::string( delim+1 );

      // see if there's another replacement rule following
      nextRule = strchr( delim, ',' );
      if ( nextRule ) 
	{
	/// if there is, remove it from the replacement string
	const int cplenNext = nextRule - delim - 1;
	replacement = replacement.substr( 0, cplenNext );
	nextRule++;
	} 
      else
	{
	nextRule = NULL;
	}
      
      // check for beginning-of-line token
      bool checkPrefixOnly = false;
      if ( pattern[0] == '^' ) 
	{
	checkPrefixOnly = true;
	}
      
      if ( checkPrefixOnly ) 
	{
	// Check if rule applies to given path.
	if ( path.substr( 0, pattern.length() - 1 ) == pattern.substr( 1 ) ) 
	  {
	  // Yes, it does: Substitute prefix accordingly and return pointer
	  // to buffer containing modified path.
	  buffer = buffer.replace( 0, pattern.length() - 1, replacement );
	  }
	} 
      else
	{
	// Substitute non-prefix occurences as well
	size_t found = buffer.find( pattern );
	while ( found != std::string::npos )
	  {
	  buffer = buffer.replace( found, pattern.length(), replacement );
	  found = buffer.find( pattern, found + replacement.length() ); // search after replaced string to avoid infinite recursive replacement
	  }
	}
      }
    }
  
  return buffer;
}

} // namespace cmtk
