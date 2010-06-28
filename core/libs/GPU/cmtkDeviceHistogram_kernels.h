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

#ifndef __cmtkDeviceHistogram_kernels_h_included_
#define __cmtkDeviceHistogram_kernels_h_included_

/** \addtogroup GPU */
//@{

/// Populate histogram from data, entirely on device.
void cmtkDeviceHistogramPopulate( float* histPtr, const float* dataPtr, const int* maskPtr, const float rangeFrom, const float rangeTo, int numberOfBins, int numberOfSamples );

/// Compute entropy from histogram on device.
void cmtkDeviceHistogramEntropy( float* result, const float* dataPtr, int numberOfBins );

//@}

#endif // #ifndef __cmtkDeviceHistogram_kernels_h_included_
