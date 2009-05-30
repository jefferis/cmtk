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

#ifndef __cmtkMatchedLandmarkList_h_included_
#define __cmtkMatchedLandmarkList_h_included_

#include <cmtkconfig.h>

#include <cmtkSmartPtr.h>
#include <cmtkLandmarkList.h>
#include <cmtkMatchedLandmark.h>

#include <list>

#ifdef CMTK_HAVE_VTK
#  include <vtkPoints.h>
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// List of matched landmarks in 3-D.
class MatchedLandmarkList :
  /// Inherit STL list container.
  public std::list<MatchedLandmark::SmartPtr>
{
public:
  /// This class.
  typedef MatchedLandmarkList Self;

  /// Smart pointer to MatchedLandmarkList.
  typedef SmartPointer<Self> SmartPtr;

  /// Default constructor.
  MatchedLandmarkList() {}

  /// Initialize from two separate landmark lists.
  MatchedLandmarkList( const LandmarkList* sourceList, const LandmarkList* targetList )
  {
    this->AddLandmarkLists( sourceList, targetList );
  }
  
  /// Initialize from two separate landmark lists.
  void  AddLandmarkLists( const LandmarkList* sourceList, const LandmarkList* targetList );
  
  /// Find landmark by name.
  SmartPointer<MatchedLandmark> FindByName( const char* name );
  
  /// Find landmark by name and return constant pointer.
  const SmartPointer<MatchedLandmark> FindByName( const char* name ) const;
  
#ifdef CMTK_HAVE_VTK
  /// Get landmarks source locations as vtkPoints list.
  vtkPoints* GetVtkPointsSource() const;

  /// Get landmarks source locations as vtkPoints list.
  vtkPoints* GetVtkPointsTarget() const;

  /** Get matched landmarks as two vtkPoints lists.
   */
  void GetMatchedVtkPoints( vtkPoints*& sourcePoints, vtkPoints*& targetPoints ) const;
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMatchedLandmarkList_h_included_
