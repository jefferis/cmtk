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

#ifndef __cmtkClassStreamStudyList_h_included_
#define __cmtkClassStreamStudyList_h_included_

#include <cmtkconfig.h>

#include <IO/cmtkStudyList.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** StudyList with class stream file system interface.
 */
class ClassStreamStudyList 
{
public:
  /// Read studylist.
  static StudyList* Read( const char *path );

  /// Merge with duplicated studylist.
  static StudyList* Merge( const StudyList& studyList, const char *path );
  /// Merge studylist in place.
  static StudyList* Merge( StudyList* const studyList, const char *path );

  /// Write studylist.
  static void Write( const char *path, const StudyList *studyList );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkClassStreamStudyList_h_included_
