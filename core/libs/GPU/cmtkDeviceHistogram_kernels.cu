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

#include <cstdio>
#include <cstdlib>

__global__
void 
cmtkDeviceHistogramEntropyKernel( float* result, const float *dataPtr )
{
  int tx = threadIdx.x;

  // first, load data into shared memory
  __shared__ float working[512]; // allocate maximum we possibly need
  working[tx] = dataPtr[tx];
  __syncthreads();

  // second, compute sum of all bin values via butterfly
  for ( int bit = 1; bit < blockDim.x; bit <<= 1 )
    {
      const float sum = working[tx] + working[tx^bit];
      __syncthreads();
      working[tx] = sum;
      __syncthreads();
    }
  
  // third, normalize
  if ( working[tx] )
    {
      working[tx] = dataPtr[tx] / working[tx];
    }
  
  // fourth, do p*log(p)
  if ( working[tx] > 0 )
    {
      working[tx] *= log( working[tx] );
    }
  else
    {
      working[tx] = 0;
    }
  __syncthreads();

  // fifth, another butterfly to compute \sum[p*log(p)]
  for ( int bit = 1; bit < blockDim.x; bit <<= 1 )
    {
      const float sum = working[tx] + working[tx^bit];
      __syncthreads();
      working[tx] = sum;
      __syncthreads();
    }

  result[tx] = -working[tx];
}

void
cmtkDeviceHistogramEntropy( float* result, const float* dataPtr, int numberOfBins )
{
  dim3 dimBlock( numberOfBins, 1 );
  dim3 dimGrid( 1, 1 );
  
  cmtkDeviceHistogramEntropyKernel<<<dimGrid,dimBlock>>>( result, dataPtr );
}

__global__
void 
cmtkDeviceHistogramPopulateKernel( float* histPtr, const float *dataPtr, const int numberOfBins )
{
  int tx = threadIdx.x;
}

void
cmtkDeviceHistogramPopulate( float* histPtr, const float* dataPtr, int numberOfBins, int numberOfSamples )
{
  // first, clear histogram memory on device
  if ( cudaMemset( histPtr, 0, sizeof( float ) * numberOfBins ) != cudaSuccess )
    {
      fputs( "ERROR: cudaMemset() failed.\n", stderr );
      exit( 1 );
    }

  // how many local copies of the histogram can we fit in shared memory?
  int device;
  cudaDeviceProp dprop;
  if ( (cudaGetDevice( &device ) != cudaSuccess) || (cudaGetDeviceProperties( &dprop, device ) != cudaSuccess ) )
    {
      fputs( "ERROR: could not get device properties.\n", stderr );
      exit( 1 );
    }
  
  int nLocalHistogramsSharedMemory = dprop.sharedMemPerBlock / (sizeof(float) * numberOfBins);
  if ( nLocalHistogramsSharedMemory > 512 )
    nLocalHistogramsSharedMemory = 512;

  dim3 dimBlock( nLocalHistogramsSharedMemory, 1 );
  dim3 dimGrid( numberOfSamples / nLocalHistogramsSharedMemory, 1 );
  
  cmtkDeviceHistogramPopulateKernel<<<dimGrid,dimBlock>>>( histPtr, dataPtr, numberOfBins );  

  const int residualSamples = numberOfSamples - nLocalHistogramsSharedMemory * (numberOfSamples / nLocalHistogramsSharedMemory);
  if ( residualSamples )
    {
      dim3 dimBlock( residualSamples, 1 );
      dim3 dimGrid( 1, 1 );
      
      cmtkDeviceHistogramPopulateKernel<<<dimGrid,dimBlock>>>( histPtr, dataPtr + numberOfSamples - residualSamples, numberOfBins );
    }
}
