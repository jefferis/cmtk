/*
//
//  Copyright 2010-2011 SRI International
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

#include "cmtkMetaInformationObject.h"

const std::string& 
cmtk::MetaInformationObject::GetMetaInfo( const std::string& key, const std::string& defaultVal ) const
{
  Self::KeyValueMapType::const_iterator it = this->m_MetaInformation.find( key );
  if ( it != this->m_MetaInformation.end() )
    return it->second;
  else
    return defaultVal;
}

void
cmtk::MetaInformationObject::SetMetaInfo( const std::string& key, const std::string& value )
{
  this->m_MetaInformation[key] = value;
}
