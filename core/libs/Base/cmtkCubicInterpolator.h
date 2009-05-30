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

#ifndef __cmtkCubicInterpolator_h_included_
#define __cmtkCubicInterpolator_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{
namespace 
Interpolators
{

/// Cubic interpolator.
class Cubic
{
public:
  /// Size of the interpolation region in grid points to the left and right.
  static const int RegionSizeLeftRight = 2;

  /// Get specific interpolation weight for relative coordinate.
  static Types::Coordinate GetWeight( const int weight, const Types::Coordinate x)
  {
    const Types::Coordinate xsquare = x * x;
    const Types::Coordinate xcube = xsquare * x;
    switch (weight)
      {
      case -1:
	return -0.5 * xcube + xsquare - 0.5 * x;
      case 0:
	return 1.5 * xcube - 2.5 * xsquare + 1;
      case 1:
	return -1.5 * xcube + 2 * xsquare + 0.5 * x;
      case 2:
	return 0.5 * xcube - 0.5 * xsquare;
      default:
#ifdef DEBUG
	std::cerr << "weight=" << weight << " shouldn't happen!" << std::endl;
	exit( 1 );
#endif
	break;
      }
    return 0;
  }
};

} // namespace Interpolators

//@}

} // namespace cmtk

#endif
