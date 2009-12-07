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

#ifndef __cmtkInterpolator_h_included_
#define __cmtkInterpolator_h_included_

#include <cmtkconfig.h>

namespace cmtk
{

/** \addtogroup Base */
//@{
/// Interpolation kernels to be used with UniformVolumeInterpolator.
namespace Interpolators
{
/// Constants for interpolation modes.
typedef enum 
{
  /// Nearest neighbour interpolation.
  NEAREST_NEIGHBOR,
  /// Partial volume interpolation.
  PARTIALVOLUME,
  /// (Tri-)linear interpolation.
  LINEAR,
  /// (Tri-)cubic interpolation.
  CUBIC,
  /// Sinc interpolation with cosine window.
  COSINE_SINC,
  /// Sinc interpolation with Hamming window.
  HAMMING_SINC,
  /// Default/unknown interpolation.
  DEFAULT
} InterpolationEnum;

} // namespace Interpolators
} // namespace cmtk

#endif // #ifndef __cmtkInterpolator_h_included_

