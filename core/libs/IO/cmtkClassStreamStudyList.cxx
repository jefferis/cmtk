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

#include "cmtkClassStreamStudyList.h"

#include <System/cmtkMountPoints.h>

#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkClassStreamAffineXform.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

void 
ClassStreamStudyList::Write
( const char *path, const StudyList* studyList )
{
  ClassStreamOutput stream;

  stream.Open( path, "studylist", ClassStreamOutput::MODE_WRITE );
  if ( stream.IsValid() ) 
    {
    StudyList::const_iterator it = studyList->begin();
    while ( it != studyList->end() ) 
      {
      stream.Begin( "source" );
      stream.WriteString( "studyname", it->first->GetFileSystemPath() );
      stream.End();
      ++it;
      }
    stream.Close();
    }
  else
    {
    StdErr << "ERROR: could not open archive " << path << "/studylist\n";
    }  
  
  stream.Open( path, "registration", ClassStreamOutput::MODE_WRITE_ZLIB );
  if ( stream.IsValid() ) 
    {
    StudyList::const_iterator it = studyList->begin();
    while ( it != studyList->end() ) 
      {
      StudyToXform targetList = it->second;
      
      std::map<Study::SmartPtr,bool> seen;
      
      StudyToXform::const_iterator tit;
      for ( tit = targetList.begin(); tit != targetList.end(); ++tit ) 
	{
	if ( seen.find( tit->first ) == seen.end() )
	  {
	  seen[tit->first] = true;
	  
	  stream.Begin( "registration" );
	  stream.WriteString( "reference_study", it->first->GetFileSystemPath() );
	  stream.WriteString( "floating_study", tit->first->GetFileSystemPath() );
	  
	  StudyToXform::const_iterator tit2;
	  for ( tit2 = targetList.begin(); tit2 != targetList.end(); ++tit2 ) 
	    {
	    if ( tit2->first == tit->first )
	      {
	      Xform::SmartPtr xform = tit2->second;
	      
	      AffineXform::SmartPtr affine = AffineXform::SmartPtr::DynamicCastFrom( xform );
	      if ( affine )
		{
		stream << (*affine);
		}
	      
	      WarpXform::SmartPtr warp = WarpXform::SmartPtr::DynamicCastFrom( xform );
	      if ( warp ) 
		{
		stream << warp;
		}
	      }
	    }
	  stream.End();
	  }
	}
      ++it;
      }
    
    stream.Close();
    }
  else
    {
    StdErr << "ERROR: could not open archive " << path << "/registration\n";
    }  
}

} // namespace cmtk
