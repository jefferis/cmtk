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

#ifndef __cmtkLandmarkPair_h_included_
#define __cmtkLandmarkPair_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkLandmark.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Matched landmark (landmark with source and target location).
class LandmarkPair :
  /// Inherit single landmark.
  public Landmark
{
public:
  /// This class.
  typedef LandmarkPair Self;

  /// Smart pointer to LandmarkPair.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const LandmarkPair.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Constructor.
  LandmarkPair( const Landmark& landmark, const Self::SpaceVectorType& target )
    : Landmark( landmark ),
      m_TargetLocation( target )
  {}

  /// Constructor.
  LandmarkPair( const std::string& name, const Self::SpaceVectorType& source, const Self::SpaceVectorType& target )
    : Landmark( name, source ),
      m_TargetLocation( target )
  {}

  /// Coordinates of this landmark.
  Self::SpaceVectorType m_TargetLocation;
};

/// Stream output operator.
std::ostream& operator<<( std::ostream& stream, const LandmarkPair& pair )
{
  stream << pair.m_Location << "\t" << pair.m_TargetLocation << "\t" << pair.m_Name << std::endl;
  return stream;
}

/// Stream input operator.
std::istream& operator>>( std::istream& stream, LandmarkPair& pair )
{
  stream >> pair.m_Location >> pair.m_TargetLocation >> pair.m_Name;
  return stream;
}

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLandmarkPair_h_included_
