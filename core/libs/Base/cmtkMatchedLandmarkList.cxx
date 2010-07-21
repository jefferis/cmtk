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

#include "cmtkMatchedLandmarkList.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
MatchedLandmarkList::AddLandmarkLists
( const LandmarkList* sourceList, const LandmarkList* targetList )
{
  LandmarkList::const_iterator it = sourceList->begin();
  while ( it != sourceList->end() ) 
    {
    const Landmark* targetLM = targetList->FindByName( (*it)->GetName() ).GetConstPtr();
    if ( targetLM ) 
      {
      SmartPointer<MatchedLandmark> newMatchedLM( new MatchedLandmark );
      newMatchedLM->SetName( (*it)->GetName() );      
      newMatchedLM->SetLocation( (*it)->GetLocation() );
      newMatchedLM->SetTargetLocation( targetLM->GetLocation() );
      this->push_back( newMatchedLM );
      }
    ++it;
    }
}

const SmartPointer<MatchedLandmark> 
MatchedLandmarkList::FindByName( const char* name ) const
{
  const_iterator it = this->begin();
  while ( it != this->end() ) 
    {
    if ( !strcmp( (*it)->GetName(), name ) ) 
      {
      return *it;
      }
    ++it;
    }
  
  return SmartPointer<MatchedLandmark>( NULL );
}

SmartPointer<MatchedLandmark> 
MatchedLandmarkList::FindByName( const char* name )
{
  iterator it = this->begin();
  while ( it != this->end() ) 
    {
    if ( !strcmp( (*it)->GetName(), name ) ) 
      {
      return *it;
      }
    ++it;
    }
  
  return SmartPointer<MatchedLandmark>( NULL );
}

} // namespace cmtk
