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

#include <cmtkStudyList.h>

#include <cmtkConsole.h>

#include <cmtkSplineWarpXform.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

const Study*
StudyList::GetStudy( const unsigned int studyIndex ) const
{
  if ( studyIndex < this->size() ) {
    const_iterator it = this->begin();
    for ( unsigned int i = 0; i < studyIndex; ++i ) ++it;
    return it->first;
  } else
    return NULL;
}

Study::SmartPtr
StudyList::GetStudy( const unsigned int studyIndex )
{
  if ( studyIndex < this->size() ) {
    const_iterator it = this->begin();
    for ( unsigned int i = 0; i < studyIndex; ++i ) ++it;
    return it->first;
  } else
    return Study::SmartPtr::Null;
}

const Study*
StudyList::FindStudyPath( const char *fileSystemPath ) const
{
  if ( ! fileSystemPath ) return NULL;

  const_iterator it = this->begin();
  while ( it != this->end() ) {
    if ( ! strcmp( it->first->GetFileSystemPath(), fileSystemPath ) )
      return it->first;
    ++it;
  }
  
  // not found: return NULL;
  return NULL;
}

Study::SmartPtr
StudyList::FindStudyPath( const char *fileSystemPath, const bool create )
{
  if ( ! fileSystemPath ) return Study::SmartPtr::Null;

  iterator it = this->begin();
  while ( it != this->end() ) {
    if ( ! strcmp( it->first->GetFileSystemPath(), fileSystemPath ) )
      return it->first;
    ++it;
  }
  
  // not found: return NULL or create;
  if ( !create )
    return Study::SmartPtr::Null;
  
  Study::SmartPtr newStudy;
  newStudy->SetFileSystemPath( fileSystemPath );
  this->AddStudy( newStudy );
  return newStudy;
}

const Study*
StudyList::FindStudyName( const char *name ) const
{
  if ( ! name ) return NULL;

  const_iterator it = this->begin();
  while ( it != this->end() ) {
    if ( ! strcmp( it->first->GetName(), name ) )
      return it->first;
    ++it;
  }
  
  // not found: return NULL;
  return NULL;
}

Study::SmartPtr
StudyList::FindStudyName( const char *name )
{
  if ( ! name ) return Study::SmartPtr::Null;

  iterator it = this->begin();
  while ( it != this->end() ) {
    if ( ! strcmp( it->first->GetName(), name ) )
      return it->first;
    ++it;
  }
  
  // not found: return NULL;
  return Study::SmartPtr::Null;
}

Study::SmartPtr
StudyList::AddStudy( const char *fileSystemPath )
{
  if ( ! fileSystemPath ) return Study::SmartPtr::Null;

  const_iterator it = this->begin();
  while ( it != this->end() ) 
    {
    // if this study is already in the list, we're done.
    if ( ! strcmp( it->first->GetFileSystemPath(), fileSystemPath ) )
      return Study::SmartPtr::Null;
    ++it;
    }
  
  Study::SmartPtr newStudy( Study::Read( fileSystemPath ) );
  if ( newStudy ) 
    {
    int suffix = 0;
    while ( this->FindStudyName( newStudy->GetName() ) ) {
    newStudy->SetMakeName( NULL, suffix++ );
    }
    
    (*this)[newStudy];
  }

  return newStudy;
}

void 
StudyList::AddStudy( Study::SmartPtr& study )
{
  if ( study.IsNull() ) return;

  const char* newStudyPath = study->GetFileSystemPath();

  const_iterator it = this->begin();
  while ( it != this->end() ) {
    // if this study is already in the list, we're done.
    if ( ! strcmp( it->first->GetFileSystemPath(), newStudyPath ) )
      return;
    ++it;
  }
  
  // insert new study into map.
  (*this)[study];
}

void 
StudyList::AddXform
( const char *fromStudyPath, const char *toStudyPath, AffineXform::SmartPtr& affineXform, WarpXform::SmartPtr& warpXform )
{
  Study::SmartPtr fromStudy = this->FindStudyPath( fromStudyPath, true /*create*/ );
  Study::SmartPtr toStudy = this->FindStudyPath( toStudyPath, true /*create*/ );

  this->AddXform( fromStudy, toStudy, affineXform, warpXform );
}

void 
StudyList::AddXform
( Study::SmartPtr& fromStudy, Study::SmartPtr& toStudy, AffineXform::SmartPtr& affineXform, WarpXform::SmartPtr& warpXform )
{
  if ( fromStudy.IsNull() || toStudy.IsNull() ) return;

  if ( affineXform ) 
    {
    Xform::SmartPtr xform = affineXform;
    (*this)[fromStudy].insert( std::multimap<Study::SmartPtr,Xform::SmartPtr>::value_type( toStudy, xform ) );
    }
  if ( warpXform ) 
    {
    Xform::SmartPtr xform = warpXform;
    (*this)[fromStudy].insert( std::multimap<Study::SmartPtr,Xform::SmartPtr>::value_type( toStudy, xform ) );
    }
}

void 
StudyList::DeleteStudy( const Study* study )
{
  iterator it = this->begin();
  while ( it != this->end() ) {
    if ( it->first == study ) {
      this->erase( it );
    }
    break;
  }
}

void
StudyList::Print() const
{
  const_iterator it = this->begin();
  while ( it != this->end() ) 
    {
    StdErr.printf( "Reference study: %s\n", it->first->GetFileSystemPath() );

    StudyToXform::const_iterator it2 = it->second.begin();
    while ( it2 != it->second.end() ) 
      {
      StdErr.printf( "\tFloating study: %s\n", it2->first->GetFileSystemPath() );

      if ( AffineXform::SmartPtr::DynamicCastFrom( it2->second ) )
	StdErr << "\t\tAffine transformation\n";
      
      if ( SplineWarpXform::SmartPtr::DynamicCastFrom( it2->second ) )
	StdErr << "\t\tSpline warp transformation\n";
      
      ++it2;
      }
    
    ++it;
    }
}

std::list<Study::SmartPtr>
StudyList::GetTargetStudyList( const Study* referenceStudy )
{
  std::list<Study::SmartPtr> result;
  
  for ( StudyList::iterator it = this->begin(); it != this->end(); ++it ) {
    if ( it->first.GetPtr() == referenceStudy ) {
      StudyToXform& stx = it->second;
      
      for ( StudyToXform::iterator xit = stx.begin(); xit != stx.end();
	      ++xit ) {
	result.push_back( xit->first );
      }
      break;
    }
  }
  return result;
}

std::list<Study::SmartPtr>
StudyList::GetIndependentStudyList( const Study* referenceStudy )
{
  std::list<Study::SmartPtr> result;
  
  for ( StudyList::iterator it = this->begin(); it != this->end(); ++it ) {
    if ( it->first.GetPtr() == referenceStudy ) {
      StudyToXform& stx = it->second;
      
      for ( StudyList::iterator it2 = this->begin(); it2 != this->end(); 
	    ++it2 ) {
	if ( stx.find( it2->first) == stx.end() )
	  result.push_back( it2->first );
      }
      break;
    }
  }
  return result;
}

bool
StudyList::Connected
( const Study::SmartPtr& s1, const Study::SmartPtr& s2 )
{
  StudyList::iterator it = this->find( s1 );

  if ( it == this->end() ) return false;

  StudyToXform::iterator xit = it->second.find( s2 );
  return (xit != it->second.end() );
}

} // namespace cmtk
