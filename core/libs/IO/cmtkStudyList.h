/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#ifndef __cmtkStudyList_h_included_
#define __cmtkStudyList_h_included_

#include <cmtkconfig.h>

#include <IO/cmtkStudy.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>

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
  const Study *FindStudyPath( const std::string& fileSystemPath ) const;

  /** Find non-constant Study object by file system path.
    */
  Study::SmartPtr FindStudyPath( const std::string& fileSystemPath /*!< Path of study to find in filesystem */, const bool create = false /*!< Flag whether to create a study that does not exist already */ );

  /// Find constant Study object by file system path.
  const Study *FindStudyName( const std::string& name /*!< Name of the study to find.*/ ) const;

  /// Find non-constant Study object by file system path.
  Study::SmartPtr FindStudyName( const std::string& name /*!< Name of the study to find.*/ );

  /// Add a new study entry.
  Study::SmartPtr AddStudy( const std::string& fileSystemPath /*!< Name of the study to add.*/ );
  
  /// Add an existing study object.
  void AddStudy( Study::SmartPtr& study /*!< Existing study object to add to the study list.*/ );
  
  /// Add a coordinate transformation between two studies.
  void AddXform( Study::SmartPtr& fromStudy, Study::SmartPtr& toStudy, AffineXform::SmartPtr& affineXform, WarpXform::SmartPtr& warpXform = WarpXform::SmartPtr::Null() );

  /// Add a coordinate transformation between two studies.
  void AddXform( const std::string& fromStudyPath /*!< Path of the study that the transformation maps from.*/, 
		 const std::string& toStudyPath /*!< Path of the study that the transformation maps to.*/, 
		 AffineXform::SmartPtr& affineXform /*!< Affine coordinate transformation between the given studies. */, 
		 WarpXform::SmartPtr& warpXform = WarpXform::SmartPtr::Null() /*!< Optional nonrigid transformation between the two given studies. */);
  
  /// Remove and delete given study object.
  void DeleteStudy( const Study* study /*!< Study object to delete from the list.*/ );

  /// Remove and delete given study object.
  void DeleteStudy( const unsigned int studyIndex /*!< Index of study to delete from the list.*/ );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStudyList_h_included_
