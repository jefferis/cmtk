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

#ifndef __cmtkDeviceImageConvolution_kernels_h_included_
#define __cmtkDeviceImageConvolution_kernels_h_included_

#include <cmtkconfig.h>

/** \addtogroup GPU */
//@{

/** Convolution of a 3D image (CUDA array) with a separable 3D kernel.
 *\warning Even though the convolution result is ultimately stored out-of-place in the
 * given target memory, the input array's content is destroyed in the process.
 */
void
cmtkDeviceImageConvolution( float* dest, const int* dims3, void* array, const int kernelLengthX, const float* kernelX, const int kernelLengthY, const float* kernelY, const int kernelLengthZ, const float* kernelZ );

/** In-place convolution of a 3D image (CUDA array) with a separable 3D kernel.
 */
void
cmtkDeviceImageConvolutionInPlace( const int* dims3, void* array, const int kernelLengthX, const float* kernelX, const int kernelLengthY, const float* kernelY, const int kernelLengthZ, const float* kernelZ );

//@}

#endif // #ifndef __cmtkDeviceImageConvolution_kernels_h_included_
