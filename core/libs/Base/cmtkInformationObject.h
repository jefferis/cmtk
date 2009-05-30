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

#ifndef __cmtkInformationObject_h_included_
#define __cmtkInformationObject_h_included_

#include <map>
#include <string>

#define CMTK_META_FS_PATH "FILESYSTEM_PATH"
#define CMTK_META_FILEFORMAT_ORIGINAL "FILEFORMAT_ORIGINAL"

#define CMTK_META_SPACE "SPACE"
#define CMTK_META_SPACE_ORIGINAL "SPACE_ORIGINAL"
#define CMTK_META_SPACE_UNITS_STRING "SPACE_UNITS_STRING"
#define CMTK_META_EXTERNAL_SPACE_ID "SPACE_ID_EXTERNAL"

#define CMTK_META_IMAGE_ORIENTATION "IMAGE_ORIENTATION"
#define CMTK_META_IMAGE_ORIENTATION_ORIGINAL "IMAGE_ORIENTATION_ORIGINAL"

#define CMTK_META_IMAGE_DIRECTION_VECTORS "IMAGE_DIRECTION_VECTORS"

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Meta-information associated with library objects.
class InformationObject
{
public:
  /// Virtual destructor template.
  virtual ~InformationObject() {};

  /// Check whether a key exists.
  bool MetaKeyExists( const std::string key ) const
  {
    return this->m_MetaInformation.find( key ) != this->m_MetaInformation.end();
  }

  /// The actual table of meta data: maps keys to values.
  mutable std::map<std::string,std::string> m_MetaInformation;
};

//@}

} // namespace cmtk

#endif
