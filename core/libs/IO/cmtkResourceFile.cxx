/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkResourceFile.h>

#include <string.h>
#include <limits.h>
#include <fstream>
#include <algorithm>

#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

void
ResourceFile::Read( const char* fileName )
{
  std::ifstream stream( fileName );
  if ( ! stream ) 
    {
    StdErr.printf( "Could not open resource file %s for reading.", fileName );
    return;
    }
  
  char line[PATH_MAX];
  std::string section = "";
  
  while ( !stream.eof() ) 
    {
    stream.getline( line, PATH_MAX );
    
    if ( ! line[0] ) continue;
    
    char* lastNonWS = line + strlen( line ) - 1;
    while ( (line <= lastNonWS) && isspace( *lastNonWS ) )
      --lastNonWS;
    *(lastNonWS+1) = 0;
    
    if ( (line[0] == '[') && (*lastNonWS == ']') ) 
      {
      *lastNonWS = 0;
      section = line+1;
      } 
    else
      {
      if ( section.length() ) 
	{
	(*this)[section].push_back( line );
	}
      }
    
    line[0] = 0; // this fixes a bug in fgets() with empty lines
    }
}

void
ResourceFile::Write( const char* fileName ) const
{
  std::ofstream stream( fileName );
  if ( ! stream ) 
    {
    StdErr.printf( "Could not open resource file %s for writing.", fileName );
    return;
    }
  
  const_iterator sectionIt = this->begin();
  while ( sectionIt != this->end() ) 
    {
    stream << "[" << sectionIt->first << "]\n";
    
    ResourceSection::const_iterator strIt = sectionIt->second.begin();
    while ( strIt != sectionIt->second.end() ) 
      {
      stream << *strIt << "\n";	
      ++strIt;
      }
    
    ++sectionIt;
    }
}

unsigned int
ResourceFile::AddUnique( const char* section, const char* entry, const unsigned int maxItems )
{
  ResourceSection& sec = (*this)[section];
  
  // Find (and delete) duplicates.
  ResourceSection::iterator it;
  while ( (it = std::find( sec.begin(), sec.end(), entry )) != sec.end() ) 
    {
    sec.erase( it );
    }
  
  // Add new entry to front of list.
  sec.push_front( entry );
  
  // If maximum length is given, trim the list.
  if ( maxItems ) 
    {
    if ( sec.size() > maxItems ) 
      {
      ResourceSection::iterator it = sec.begin();
      for ( unsigned int i = 0; i < maxItems; ++i ) ++it;
      sec.erase( it, sec.end() );
      }
    }
  return sec.size();
}

} // namespace cmtk
