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

#include "cmtkClassStream.h"

namespace
cmtk
{

/** \addtogroup IO */
//@{

ClassStream& 
ClassStream::operator << ( const LandmarkList *landmarkList )
{
  if ( ! landmarkList ) return *this;

  for ( LandmarkList::const_iterator it = landmarkList->begin(); it != landmarkList->end(); ++it )
    {
    this->Begin( "landmark" );
    this->WriteString( "name", (*it)->m_Name.c_str() );
    this->WriteCoordinateArray( "location", (*it)->m_Location.begin(), (*it)->m_Location.Size() );
    this->End();
  }

  return *this;
}

ClassStream& 
ClassStream::operator >> ( LandmarkList::SmartPtr& landmarkList )
{
  if ( !this->IsValid() ) 
    {
    landmarkList = LandmarkList::SmartPtr::Null();
    return *this;
    }
  
  landmarkList = LandmarkList::SmartPtr( new LandmarkList );
  
  while ( this->Seek( "landmark" ) ) 
    {
    Landmark::SmartPtr landmark( new Landmark );
    
    char *tmpStr = this->ReadString( "name" );
    if ( tmpStr ) 
      {
      landmark->m_Name = tmpStr;
      free( tmpStr );
      }
    
    Types::Coordinate location[3];
    this->ReadCoordinateArray( "location", location, 3 );
    landmark->m_Location = Landmark::SpaceVectorType( location );
    
    landmarkList->insert( landmarkList->end(), landmark );
    }
  
  return *this;
}

} // namespace cmtk
