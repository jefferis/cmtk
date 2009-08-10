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

#ifndef __cmtkLandmarkListVTK_h_included_
#define __cmtkLandmarkListVTK_h_included_

#include <cmtkconfig.h>

#include <cmtkSmartPtr.h>
#include <cmtkLandmarkList.h>

#include <list>

#include <vtkPoints.h>

namespace
cmtk
{

/** \addtogroup VTKWrapper */
//@{

/// List of landmarks in 3-D with VTK support.
class LandmarkListVTK :
  /// Inherit non-VTK base class
  public LandmarkList
{
public:
  /// This class.
  typedef LandmarkListVTK Self;

  /// Smart pointer to LandmarkListVTK.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef LandmarkList Superclass;

  /// Copy-and-convert constructor from superclass.
  LandmarkListVTK( const Superclass& other ) : Superclass( other ) {};

  /// Get landmarks as vtkPoints list.
  vtkPoints* GetVtkPoints() const;

  /** Get matched landmarks as two vtkPoints lists.
   * The landmarks in this list and the list given as the "targetLL" parameter
   * are matched using their name and inserted into two vtkPoints objects. The
   * order of insertion is determined by their order in this object. Landmarks
   * that have no equally named counterpart in the other list will be ignored.
   */
  vtkPoints* GetMatchedVtkPoints( vtkPoints*& targetPoints, const LandmarkListVTK *targetLL ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLandmarkListVTK_h_included_
