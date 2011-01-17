/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkTypedStreamStudylist_h_included_
#define __cmtkTypedStreamStudylist_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>

#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Studylist with typedstream file system interface.
 * This class provides the necessary functions to read studylist objects
 * from typedstream archives.
 */
class TypedStreamStudylist 
{
public:
  /// Smart pointer to TypedStreamStudylist
  typedef SmartPointer<TypedStreamStudylist> SmartPtr;

  /// Default conctructor.
  TypedStreamStudylist();

  /// Destructor.
  ~TypedStreamStudylist();

  /** Read constructor.
   */
  TypedStreamStudylist ( const char *studylistpath /*!<  The typedstream archive to read the object from. */) 
  { 
    this->Clear();
    this->Read( studylistpath );
  }

  /// Read object from disk.
  bool Read( const char* studylistpath );

  /// Return affine transformation as stored in the studylist.
  AffineXform::SmartPtr& GetAffineXform() 
  { 
    return this->m_AffineXform; 
  }

  /// Return local deformation as stored in the studylist.
  WarpXform::SmartPtr& GetWarpXform() 
  { 
    return this->m_WarpXform;
  }

  /// Return study path.
  const char* GetStudyPath( const int index ) const 
  {
    return StudyPath[index];
  }

  /// Return reference study path.
  const char* GetReferenceStudyPath() const 
  {
    return StudyPath[ReferenceStudyIndex];
  }

  /// Return floating study path.
  const char* GetFloatingStudyPath( const int floatingIndex = 0 ) const 
  {
    if ( floatingIndex < ReferenceStudyIndex )
      return StudyPath[floatingIndex];
    else
      if ( floatingIndex < 1 )
	return StudyPath[floatingIndex+1];
      else
	return NULL;
  }
  
  /// Create studylist archive.
  static void Write( const char *path, const char* referenceStudy, const char* floatingStudy, const Xform* xform );

private:
  /// The names of the two studies referenced in the studylist.
  char *StudyPath[2];
  
  /// Index of the reference study among the studies in this list (from 0).
  int ReferenceStudyIndex;

  /// Pointer to the affine transformation of this studylist.
  AffineXform::SmartPtr m_AffineXform;

  /// Pointer to the local deformation of this studylist.
  WarpXform::SmartPtr m_WarpXform;

  /// Initialize all internal data structures.
  void Clear();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTypedStreamStudylist_h_included_
