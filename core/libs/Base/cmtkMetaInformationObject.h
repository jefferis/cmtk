/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkMetaInformationObject_h_included_
#define __cmtkMetaInformationObject_h_included_

#include <map>
#include <string>

namespace
cmtk
{

const char* const META_FS_PATH = "FILESYSTEM_PATH";
const char* const META_FILEFORMAT_ORIGINAL = "FILEFORMAT_ORIGINAL";

const char* const META_SPACE = "SPACE";
const char* const META_SPACE_ORIGINAL = "SPACE_ORIGINAL";
const char* const META_SPACE_UNITS_STRING = "SPACE_UNITS_STRING";
const char* const META_EXTERNAL_SPACE_ID = "SPACE_ID_EXTERNAL";

const char* const META_IMAGE_ORIENTATION = "IMAGE_ORIENTATION";
const char* const META_IMAGE_ORIENTATION_ORIGINAL = "IMAGE_ORIENTATION_ORIGINAL";

const char* const META_IMAGE_DIRECTION_VECTORS = "IMAGE_DIRECTION_VECTORS";

/** \addtogroup Base */
//@{

/// Meta-information associated with library objects.
class MetaInformationObject
{
public:
  /// Default constructor: do nothing.
  MetaInformationObject() {}

  /// Copy constructor: copy meta information when copying higher-level objects.
  MetaInformationObject( const MetaInformationObject& other )
    : m_MetaInformation( other.m_MetaInformation )
  {}

  /// Virtual destructor template.
  virtual ~MetaInformationObject() {};

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
