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

#ifndef __cmtkMatchedLandmarkListVTK_h_included_
#define __cmtkMatchedLandmarkListVTK_h_included_

#include <cmtkconfig.h>

#include <cmtkSmartPtr.h>
#include <cmtkLandmarkListVTK.h>
#include <cmtkMatchedLandmark.h>
#include <cmtkMatchedLandmarkList.h>

#include <list>

#ifdef CMTK_HAVE_VTK
#  include <vtkPoints.h>
#endif

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

/// List of matched landmarks in 3-D with support for VTK data structures.
class MatchedLandmarkListVTK :
  /// Inherit from non-VTK base class.
  public MatchedLandmarkList
{
public:
  /// This class.
  typedef MatchedLandmarkListVTK Self;

  /// Parent class.
  typedef MatchedLandmarkList Superclass;

  /// Smart pointer to MatchedLandmarkListVTK.
  typedef SmartPointer<Self> SmartPtr;

  /// Initialize from two separate landmark lists.
  MatchedLandmarkListVTK( const LandmarkListVTK* sourceList, const LandmarkListVTK* targetList ) : Superclass( sourceList, targetList ) {};
  
  /// Get landmarks source locations as vtkPoints list.
  vtkPoints* GetVtkPointsSource() const;

  /// Get landmarks source locations as vtkPoints list.
  vtkPoints* GetVtkPointsTarget() const;

  /** Get matched landmarks as two vtkPoints lists.
   */
  void GetMatchedVtkPoints( vtkPoints*& sourcePoints, vtkPoints*& targetPoints ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMatchedLandmarkListVTK_h_included_
