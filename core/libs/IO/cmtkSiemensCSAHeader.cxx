/*
//
//  Copyright 2004-2012 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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
//  $Revision: 4527 $
//
//  $LastChangedDate: 2012-10-02 10:16:03 -0700 (Tue, 02 Oct 2012) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include "cmtkSiemensCSAHeader.h"

#include <IO/cmtkFileConstHeader.h>

#include <System/cmtkConsole.h>

cmtk::SiemensCSAHeader::SiemensCSAHeader( const char* csaData, const size_t csaLength )
{
  FileConstHeader fileHeader( csaData, false /*isBigEndian*/ ); // Siemens CSA header is always little endian

  const char headerID2[] = { 'S', 'V', '1', '0' };
  const bool csa2hdr = fileHeader.CompareFieldStringN( 0, headerID2, sizeof( headerID2 ) );

  // get number of items according to either CSA1 or CSA2 format
  const size_t nTags = csa2hdr ? fileHeader.GetField<unsigned int>( 8 ) : fileHeader.GetField<unsigned int>( 0 );
    
  size_t tagOffset = csa2hdr ? 16 : 8; // start after header: length is 16 bytes for CSA2, 8 bytes for CSA1
  for ( size_t tag = 0; tag < nTags; ++tag )
    {
    // first, get tag name (up to 64 characters plus 0x0)
    char tagName[65];
    fileHeader.GetFieldString( tagOffset, tagName, 64 );

    // find number of items for this tag and make room in allocated vector
    const size_t nItems = fileHeader.GetField<unsigned int>( tagOffset + 76 );

    // create new tag object
    Self::value_type newTag( tagName, std::vector<std::string>() );
    newTag.second.resize( nItems );
    
    tagOffset += 84;
    for ( size_t item = 0; item < nItems; ++item )
      {
      const size_t itemLen = fileHeader.GetField<unsigned int>( tagOffset );

      if ( itemLen )
	{
	std::vector<char> itemStr( itemLen );
	fileHeader.GetFieldString( tagOffset+16, &(itemStr[0]), itemLen );	
	newTag.second[item] = std::string( itemStr.begin(), itemStr.end() );
	}
      
      tagOffset += 4*((itemLen+3)/4) /*move up to nearest 4-byte boundary*/ + 16 /*the 4 ints at the beginning of item, including itemLength*/;
      }

    // finally, insert tag into the map
    std::pair<Self::iterator,bool> inserted = this->insert( newTag );
    if ( !inserted.second )
      {
      StdErr << "Warning: CSA tag named '" << tagName << "' appears more than once.\n";
      }
    }
}

std::ostream& 
cmtk::operator<<( std::ostream& stream, const SiemensCSAHeader& csaHeader )
{
  for ( SiemensCSAHeader::const_iterator it = csaHeader.begin(); it != csaHeader.end(); ++it )
    {
    stream << it->first << " nitems=" << it->second.size() << "\n";
    
    for ( size_t item = 0; item < it->second.size(); ++item )
      {
      stream << "\t\"" << it->second[item] << "\" [" << it->second[item].length() << "]\n" ;
      }
    }

  return stream;
}
