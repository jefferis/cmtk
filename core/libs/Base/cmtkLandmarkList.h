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

#ifndef __cmtkLandmarkList_h_included_
#define __cmtkLandmarkList_h_included_

#include <cmtkconfig.h>

#include <System/cmtkSmartPtr.h>
#include <Base/cmtkLandmark.h>

#include <list>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// List of landmarks.
class LandmarkList :
  /// Inherit STL list container.
  public std::list< SmartPointer<Landmark> >
{
public:
  /// This class.
  typedef LandmarkList Self;

  /// Smart pointer to LandmarkList.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const LandmarkList.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Find landmark by name.
  Landmark::SmartPtr FindByName( const std::string& name );

  /// Find landmark by name and return constant pointer.
  Landmark::SmartConstPtr FindByName( const std::string& name ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLandmarkList_h_included_
