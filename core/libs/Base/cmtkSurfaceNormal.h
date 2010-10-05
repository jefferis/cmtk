/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkSurfaceNormal_h_included_
#define __cmtkSurfaceNormal_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkFixedVector.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/// Class that computes the surface normal.
class
SurfaceNormal
{
public:
  /// This class.
  typedef SurfaceNormal Self;

  /// Space vector type.
  typedef FixedVector<3,Types::Coordinate> SpaceVectorType;

  /// Constructor: takes two non-collinear vectors that span the surface.
  SurfaceNormal( const SpaceVectorType& s1, const SpaceVectorType& s2 )
  {
    this->m_Normal = FixedVectorStaticInitializer<3,Types::Coordinate>::Init( s1[1] * s2[2] - s1[2] * s2[1], s1[2] * s2[0] - s1[0] * s2[2], s1[0] * s2[1] - s1[1] * s2[0] );
  }

  /// Get the normal vector.
  const SpaceVectorType& Get() const
  {
    return this->m_Normal;
  }

private:
  /// The surface normal vector.
  SpaceVectorType m_Normal;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSurfaceNormal_h_included_
