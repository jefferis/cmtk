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

#ifndef __cmtkSimpleLevelsetDevice_kernels_h_included_
#define __cmtkSimpleLevelsetDevice_kernels_h_included_

#include <cmtkconfig.h>

/** \addtogroup GPU */
//@{

namespace
cmtk
{

/// Count inside pixels and compute the sums of inside and outside pixels.
void
SimpleLevelsetDeviceUpdateInsideOutside( float* levelset /*!< Input: current levelset */, float* volume /*!< Input: image data */, const int nPixels /*!< Input: number of pixels */,
					 float* insideSum /*!< Output: sum of data values in the "inside" region */, 
					 float* outsideSum /*!< Output: sum of data values in the "outside" region */, 
					 int* nInside /*!< Output: number of pixels in the "inside" region */ );

/// Update levelset based on mean values in inside and outside regions, then threshold the result.
void
SimpleLevelsetDeviceUpdateLevelset( float* levelset /*!< Input/output: current levelset, is updated by this function */, 
				    float* volume /*!< Input: image data */, const int nPixels /*!< Input: number of pixels */,
				    const float mInside /*!< Input: mean value of data in the "inside" region */, 
				    const float mOutside /*!< Input: mean value of data in the "outside" region */, 
				    const float ratioInOut /*!< Input: ratio of inside and outside pixel counts. */,
				    const float timeDelta /*!< Input: levelset evolution time constant */,
				    const float levelsetThreshold /*!< Input: levelset threshold */ );
} // namespace cmtk

//@}

#endif // #ifndef __cmtkSimpleLevelsetDevice_kernels_h_included_
