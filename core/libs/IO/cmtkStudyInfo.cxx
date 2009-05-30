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

#include <cmtkStudyInfo.h>

#include <string.h>

#ifdef _MSC_VER
#  include <malloc.h>
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

StudyInfo LastImageStudyInfo;

/// Geometry dependent tag.
#define CMTK_TAG_Geometry 0

/// Content dependent tag.
#define CMTK_TAG_Content 1

const char StudyInfo::DefaultTagAssignmentsTable[INFO_NUMFIELDS] =
{ CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, 
  CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, 
  CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, 
  CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, 
  CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, 
  CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content, 
  CMTK_TAG_Geometry, CMTK_TAG_Content, CMTK_TAG_Content, CMTK_TAG_Content,
  CMTK_TAG_Geometry, CMTK_TAG_Geometry, CMTK_TAG_Geometry, CMTK_TAG_Geometry, 
  CMTK_TAG_Geometry, CMTK_TAG_Geometry, CMTK_TAG_Geometry 
};

StudyInfo::StudyInfo ( const StudyInfo& geometry, const StudyInfo& content, const char* assignments )
{
  Clear();
  Merge( geometry, content, assignments );
}

StudyInfo::~StudyInfo ()
{
  for ( int i = 0; i < INFO_NUMFIELDS; ++i )
    if ( Elements[i] ) 
      {
      free(Elements[i]);
      Elements[i] = NULL;
      }
}

StudyInfo* StudyInfo::Clone() const 
{
  StudyInfo *instance = new StudyInfo;
  
  instance->SliceOffset = SliceOffset;
  instance->SliceDirection = SliceDirection;

  for ( int idx = 0; idx<INFO_NUMFIELDS; ++idx ) 
    {
    if ( Elements[idx] )
      instance->Elements[idx] = strdup( Elements[idx] );
    else
      instance->Elements[idx] = NULL;
    }
  
  return instance;
}

const char* StudyInfo::GetField ( const int idx, const int anon ) const 
{
  if ( anon && (idx < INFO_NUMANONS) )
    return StrAnonymized();
  else
    if (Elements[idx])
      return Elements[idx];
    else
      return StrNotAvailable();
}

const char* StudyInfo::ReturnFieldPtr ( const int idx ) const 
{
  return Elements[idx];
}

void StudyInfo::SetField
( const int idx, char *const str, const bool dup )
{
  if ( Elements[idx] ) 
    {
    if ( str )
      if ( dup && !strcmp( str, Elements[idx] ) ) return;
    free( Elements[idx] );
    }
  
  if ( dup && str )
    Elements[idx] = strdup(str);
  else
    Elements[idx] = str;
}

StudyInfo& StudyInfo::Merge 
( const StudyInfo& geometry, const StudyInfo& content,
  const char* assignments )
{
  const char* Assign = (assignments) ? assignments : DefaultTagAssignmentsTable;

  for ( int idx = 0; idx<INFO_NUMFIELDS; ++idx ) 
    {
    if ( Elements[idx] )
      free( Elements[idx] );
    
    if ( Assign[idx] == CMTK_TAG_Geometry )
      SetField( idx, strdup( geometry.GetField( idx, 0 ) ), false /*duplicate*/ );
    else
      SetField( idx, strdup( content.GetField( idx, 0 ) ), false /*duplicate*/);
    }
  
  SliceOffset = geometry.SliceOffset;
  SliceDirection = geometry.SliceDirection;
  
  return *this;
}

StudyInfo& StudyInfo::operator= ( const StudyInfo& other )
{
  for ( int idx = 0; idx<INFO_NUMFIELDS; ++idx ) 
    {
    SetField( idx, strdup( other.GetField( idx, StudyInfo_NoAnon ) ), false /*duplicate*/ );
    }
  return *this;
}

} // namespace cmtk
