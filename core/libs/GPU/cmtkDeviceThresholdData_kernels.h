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

#ifndef __cmtkDeviceThresholdData_kernels_h_included_
#define __cmtkDeviceThresholdData_kernels_h_included_

#include <cmtkconfig.h>

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/** Threshold data on device.
 */
void
DeviceThresholdData( float* dest, const int n, const float lowerThreshold, const float upperThreshold );

} // namespace cmtk

//@}

#endif // #ifndef __cmtkDeviceThresholdData_kernels_h_included_
