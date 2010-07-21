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

#include "cmtkStrUtility.h"

#include <cstring>
#include <cstdlib>
#include <climits>

#ifdef _MSC_VER
#  include <direct.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

int
StrCmp( const char* s1, const char* s2 )
{
  if ( s1 == NULL ) 
    {
    if ( s2 == NULL ) return 0;
    else return -1;
    } 
  else
    {
    if ( s2 == NULL ) return 1;
    else return strcmp( s1, s2 );
    }
}

void
StrReplace( char*& s1, const char* s2 )
{
  StrFree( s1 );
  if ( s2 )
    s1 = strdup( s2 );
  else
    s1 = NULL;
}

void StrFree( char *const s )
{
  free( s );
}

int StrPrefixCmp( const char *s, const char* prefix )
{
  return (0 == strncmp( s, prefix, strlen( prefix ) ));
}

static char StrBuffer[PATH_MAX];

const char*
StrDir( const char *path )
{
  const char *slash = strrchr( path, CMTK_PATH_SEPARATOR );
  if ( slash && (slash != path) ) 
    {
    int dirLen = (slash-path);
    strncpy( StrBuffer, path, dirLen );
    StrBuffer[dirLen] = 0;
    } 
  else
    {
    if ( slash )
      strcpy( StrBuffer, CMTK_PATH_SEPARATOR_STR );
    else
      strcpy( StrBuffer, path );
    }
  return StrBuffer;
}

const char*
StrFName( const char *path )
{
  const char *slash = strrchr( path, CMTK_PATH_SEPARATOR );
  if ( slash )
    {
    strcpy( StrBuffer, slash+1 );
    } 
  else
    {
    StrBuffer[0] = 0;
    }
  return StrBuffer;
}

std::string
StrReplace
( const std::string& str, const std::map<std::string,std::string>& rules, const bool multiple )
{
  std::string result = str;
  
  std::map<std::string,std::string>::const_iterator it = rules.begin();
  while ( it != rules.end() ) 
    {
    bool replaced = true;
    while ( replaced ) 
      {
      replaced = false;
      std::string::size_type pos = result.find( it->first );
      while ( pos != std::string::npos ) 
	{
	result.replace( pos, it->first.length(), it->second );
	replaced = true;
	pos = result.find( it->first );
	if ( ! multiple ) break;
	}
      
      if ( ! multiple ) break;
      }
    ++it;
    }
  return result;
}

} // namespace cmtk
