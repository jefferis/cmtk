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

#include <iostream>

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
  LandmarkPair( const Landmark& landmark, const Self::SpaceVectorType& target,
		const Types::Coordinate residual = -1 /*!< Landmark matching residual (e.g., wqith respect to linear fit) */, 
		const bool precise = true /*!< Flag whether this landmark is "precise" enough for use in registration etc. */ )
    : Landmark( landmark ),
      m_TargetLocation( target ),
      m_Residual( residual ),
      m_Precise( precise )
  {}

  /// Constructor.
  LandmarkPair( const std::string& name, const Self::SpaceVectorType& source, const Self::SpaceVectorType& target,
    		const Types::Coordinate residual = -1 /*!< Landmark matching residual (e.g., wqith respect to linear fit) */, 
		const bool precise = true /*!< Flag whether this landmark is "precise" enough for use in registration etc. */ )
    : Landmark( name, source ),
      m_TargetLocation( target ),
      m_Residual( residual ),
      m_Precise( precise )
  {}

  /// Coordinates of this landmark.
  Self::SpaceVectorType m_TargetLocation;

  /// Fitting residual (negative if unknown).
  Types::Coordinate m_Residual;

  /// Precision flag. Only landmarks with this flag set should be used for registration.
  bool m_Precise;
};

/// Stream output operator.
std::ostream& operator<<( std::ostream& stream, const LandmarkPair& pair );

/// Stream input operator.
std::istream& operator>>( std::istream& stream, LandmarkPair& pair );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkLandmarkPair_h_included_
