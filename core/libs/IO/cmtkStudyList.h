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

#ifndef __cmtkStudyList_h_included_
#define __cmtkStudyList_h_included_

#include <cmtkconfig.h>

#include <cmtkStudy.h>

#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>

#include <map>
#include <list>

namespace
cmtk
{

/** \addtogroup IO */
//@{

typedef std::multimap<Study::SmartPtr, Xform::SmartPtr> StudyToXform;
typedef std::map< Study::SmartPtr, StudyToXform > StudyToStudyToXform;

class StudyList :
  /// Also inherit container functions from STL list class.
  public StudyToStudyToXform
{
public:
  /// This class.
  typedef StudyList Self;

  /// Smart pointer to StudyList.
  typedef SmartPointer<Self> SmartPtr;

  /// Parent class.
  typedef StudyToStudyToXform Superclass;

  /** Default constructor.
   * Do nothing really, except be there.
   */
  StudyList() {};

  /** Copy constructor.
   *\todo Implement copying of transformation links between studies.
   */
  StudyList( const StudyList& slist ) : Superclass( slist ) {};

  /// Get constant Study object by index.
  const Study *GetStudy( const unsigned int studyIndex ) const;

  /// Get non-constant Study object by index.
  Study::SmartPtr GetStudy( const unsigned int studyIndex );  

  /// Find constant Study object by file system path.
  const Study *FindStudyPath( const char *fileSystemPath ) const;

  /** Find non-constant Study object by file system path.
    *\param create If true, any studies not found will be created and added
    * to this list object.
    */
  Study::SmartPtr FindStudyPath( const char *fileSystemPath, const bool create = false );

  /// Find constant Study object by file system path.
  const Study *FindStudyName( const char *name ) const;

  /// Find non-constant Study object by file system path.
  Study::SmartPtr FindStudyName( const char *name );

  /// Add a new study entry.
  Study::SmartPtr AddStudy( const char *fileSystemPath );
  
  /// Add an existing study object.
  void AddStudy( Study::SmartPtr& study );
  
  /// Add a coordinate transformation between two studies.
  void AddXform( Study::SmartPtr& fromStudy, Study::SmartPtr& toStudy, AffineXform::SmartPtr& affineXform, WarpXform::SmartPtr& warpXform = WarpXform::SmartPtr::Null );

  /// Add a coordinate transformation between two studies.
  void AddXform( const char *fromStudyPath, const char *toStudyPath, AffineXform::SmartPtr& affineXform, WarpXform::SmartPtr& warpXform = WarpXform::SmartPtr::Null );
  
  /// Remove and delete given study object.
  void DeleteStudy( const Study* study );

  /// Remove and delete given study object.
  void DeleteStudy( const unsigned int studyIndex );

  /// Print studylist information.
  void Print() const;

  /// Get all target studies for a given reference study.
  std::list<Study::SmartPtr> GetTargetStudyList( const Study* referenceStudy );

  /// Get all independent (non-target) studies for a given reference study.
  std::list<Study::SmartPtr> GetIndependentStudyList( const Study* referenceStudy );

  /// Test whether two studies are connected.
  bool Connected( const Study::SmartPtr& s1, const Study::SmartPtr& s2 );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStudyList_h_included_
