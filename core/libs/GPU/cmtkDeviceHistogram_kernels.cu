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

#include "cmtkDeviceHistogram_kernels.h"

__global__
void 
cmtkDeviceHistogramPLogPKernel( float* result, const float *dataPtr )
{
  int tx = threadIdx.x;

  // first, load data into shared memory
  __shared__ float working[512]; // allocate maximum we possibly need
  working[tx] = dataPtr[tx];
  __syncthreads();

  // second, compute sum of all bin values via butterfly
  for ( int bit = 1; bit <= blockDim.x; bit <<= 1 )
    {
      working[tx] += working[tx^bit];
      __syncthreads();
    }
  
  // third, normalize
  working[tx] = dataPtr[tx] / working[tx];
  __syncthreads();

  // fourth, do p*log(p)
  working[tx] *= log( working[tx] );
  __syncthreads();

  // fifth, another butterfly to compute \sum[p*log(p)]
  for ( int bit = 1; bit <= blockDim.x; bit <<= 1 )
    {
      working[tx] += working[tx^bit];
      __syncthreads();
    }

  result[tx] = working[tx];
}

void
cmtkDeviceHistogramPLogP( float* result, const float* dataPtr, int numberOfBins )
{
  dim3 dimBlock( numberOfBins, 1 );
  dim3 dimGrid( 1, 1 );
  
  cmtkDeviceHistogramPLogPKernel<<<dimGrid,dimBlock>>>( result, dataPtr );
}
