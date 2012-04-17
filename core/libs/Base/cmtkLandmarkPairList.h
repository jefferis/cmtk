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

#ifndef __cmtkLandmarkPairList_h_included_
#define __cmtkLandmarkPairList_h_included_

#include <cmtkconfig.h>

#include <System/cmtkSmartPtr.h>
#include <Base/cmtkLandmarkList.h>
#include <Base/cmtkLandmarkPair.h>

#include <list>
#include <string>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// List of matched landmarks in 3-D.
class LandmarkPairList :
  /// Inherit STL list container.
  public std::list<LandmarkPair>
{
public:
  /// This class.
  typedef LandmarkPairList Self;

  /// Smart pointer to LandmarkPairList.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const LandmarkPairList.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// List iterator.
  typedef std::list<LandmarkPair>::iterator Iterator;

  /// List const iterator.
  typedef std::list<LandmarkPair>::const_iterator ConstIterator;

  /// Default constructor.
  LandmarkPairList() {}

  /// Initialize from two separate landmark lists.
  LandmarkPairList( const LandmarkList& sourceList, const LandmarkList& targetList )
  {
    this->AddLandmarkLists( sourceList, targetList );
  }
  
  /// Initialize from two separate landmark lists.
  void AddLandmarkLists( const LandmarkList& sourceList, const LandmarkList& targetList );
  
  /// Find landmark by name.
  Self::Iterator FindByName( const std::string& name );
  
  /// Find landmark by name and return constant pointer.
  Self::ConstIterator FindByName( const std::string& name ) const;
};

/// Stream output operator.
std::ostream& operator<<( std::ostream& stream, const LandmarkPairList& pairList )
{
  for ( LandmarkPairList::ConstIterator it = pairList.begin(); it != pairList.end(); ++it )
    stream << it->m_Location << "\t" << it->m_TargetLocation << "\t" << it->m_Name << std::endl;
  return stream;
}

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLandmarkPairList_h_included_
