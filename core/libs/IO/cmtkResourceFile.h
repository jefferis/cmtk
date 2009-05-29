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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkResourceFile_h_included_
#define __cmtkResourceFile_h_included_

#include <cmtkconfig.h>

#include <string>
#include <map>
#include <list>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Resource file section is a list of strings.
typedef std::list<std::string> ResourceSection;

/** Resource file is a map of sections accessed by section title.
 */
class ResourceFile :
  /// Resource file is a map from section titles to resource sections.
  public std::map< std::string, ResourceSection >
{
public:
  /// Default constructor: do nothing.
  ResourceFile() {};

  /// Read constructor.
  ResourceFile( const char* fileName ) 
  {
    this->Read( fileName );
  }

  /// Read from resource file.
  void Read( const char* fileName );

  /// Write to resource file.
  void Write( const char* fileName ) const;

  /** Add a unique (non-duplicate) string to resource section.
   *\return Number of entries that are in the section after adding new entry.
   */
  unsigned int AddUnique( const char* section, const char* entry, const unsigned int maxItems = 0 );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkResourceFile_h_included_
