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

#include "cmtkStrUtility.h"

#include <string.h>

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

const char*
StrNStr( const char* haystack, const size_t nBytes, const char* needle )
{
  for ( size_t hofs = 0; hofs < nBytes; ++hofs )
    {
    size_t hidx = hofs;
    const char* nchr = needle;

    while ( (*nchr!=0) && (hidx<nBytes) && (*nchr==haystack[hidx]) ) 
      {
      ++nchr;
      ++hidx;
      }
    
    // found?
    if ( *nchr == 0 )
      return haystack+hofs;
    }
  
  return NULL;
}

int StrPrefixCmp( const char *s, const char* prefix )
{
  return (0 == strncmp( s, prefix, strlen( prefix ) ));
}

std::string
StrReplaceByRules
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

std::string
StrReplace( const std::string& str, const std::string& search, const std::string& replace )
{
  std::string result = str;

  if ( ! search.empty() )
    {
    for ( size_t b = result.find(search, 0); b != result.npos; b = result.find(search, b) )
      {
      result.replace( b, search.size(), replace );
      b += (replace.size() - search.size());
      }
    }

  return result;
}

std::string
StrMakeLegalInPath( const std::string& s )
{
  std::string result = s;
  
  result = StrReplace( result, " ", "_" );      
  result = StrReplace( result, ":", "_" );  
  
  return result;
}


} // namespace cmtk
