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

#ifndef __cmtkLandmarkList_h_included_
#define __cmtkLandmarkList_h_included_

#include <cmtkconfig.h>

#include <cmtkSmartPtr.h>
#include <cmtkLandmark.h>

#include <list>

#ifdef CMTK_HAVE_VTK
#  include <vtkPoints.h>
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// List of landmarks in 3-D.
class LandmarkList :
  /// Inherit STL list container.
  public std::list< SmartPointer<Landmark> >
{
public:
  /// This class.
  typedef LandmarkList Self;

  /// Smart pointer to LandmarkList.
  typedef SmartPointer<Self> SmartPtr;

  /// Find landmark by name.
  Landmark::SmartPtr FindByName( const char* name );

  /// Find landmark by name and return constant pointer.
  const Landmark::SmartPtr FindByName( const char* name ) const;

#ifdef CMTK_HAVE_VTK
  /// Get landmarks as vtkPoints list.
  vtkPoints* GetVtkPoints() const;

  /** Get matched landmarks as two vtkPoints lists.
   * The landmarks in this list and the list given as the "targetLL" parameter
   * are matched using their name and inserted into two vtkPoints objects. The
   * order of insertion is determined by their order in this object. Landmarks
   * that have no equally named counterpart in the other list will be ignored.
   */
  vtkPoints* GetMatchedVtkPoints( vtkPoints*& targetPoints, const LandmarkList *targetLL ) const;
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLandmarkList_h_included_
