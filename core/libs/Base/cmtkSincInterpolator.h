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

#ifndef __cmtkSincInterpolator_h_included_
#define __cmtkSincInterpolator_h_included_

#include <cmtkconfig.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{
namespace 
Interpolators
{

/// Sinc interpolator with Hamming window.
template<int NRadius=5>
class HammingSinc
{
public:
  /// This class.
  typedef HammingSinc<NRadius> Self;

  /// Size of the interpolation region in grid points to the left and right.
  static const int RegionSizeLeftRight = NRadius;

  /// Get specific interpolation weight for relative coordinate.
  static Types::Coordinate GetWeight( const int i, const Types::Coordinate x )
  {
    const Types::Coordinate piDiff = M_PI * (x - i);
#ifdef CMTK_NO_INCLASS_MEMBER_INITIALIZATION
    const Types::Coordinate result = 0.54 + 0.46 * cos( piDiff / Self::RegionSizeLeftRight ) * sin( piDiff ) / piDiff;
#else
    const Types::Coordinate result = 0.54 + 0.46 * cos( piDiff * Self::InternalFactor ) * sin( piDiff ) / piDiff;
#endif
    return finite( result ) ? result : 1;
  }

#ifndef CMTK_NO_INCLASS_MEMBER_INITIALIZATION
private:
  /// Internal factor.
  static const Types::Coordinate InternalFactor = 1.0 / Self::RegionSizeLeftRight;
#endif
};

/// Sinc interpolator with Cosine window.
template<int NRadius=5>
class CosineSinc
{
public:
  /// This class.
  typedef CosineSinc<NRadius> Self;

  /// Size of the interpolation region in grid points to the left and right.
  static const int RegionSizeLeftRight = NRadius;

  /// Get specific interpolation weight for relative coordinate.
  static Types::Coordinate GetWeight( const int i, const Types::Coordinate x )
  {
    const Types::Coordinate piDiff = M_PI * (x - i);
#ifdef CMTK_NO_INCLASS_MEMBER_INITIALIZATION
    const Types::Coordinate result = cos( piDiff / (2*Self::RegionSizeLeftRight) ) * sin( piDiff ) / piDiff;
#else
    const Types::Coordinate result = cos( piDiff * Self::InternalFactor ) * sin( piDiff ) / piDiff;
#endif
    return finite( result ) ? result : 1;
  }

#ifndef CMTK_NO_INCLASS_MEMBER_INITIALIZATION
private:
  /// Internal factor.
  static const Types::Coordinate InternalFactor = 1.0 / (2*Self::RegionSizeLeftRight);
#endif
};

} // namespace Interpolators

//@}

} // namespace cmtk

#endif
