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

extern __shared__ float shared[];

__global__
void 
cmtkDeviceHistogramEntropyKernel( float* result, const float *dataPtr )
{
  int tx = threadIdx.x;

  // first, load data into shared memory
  float *working = &shared[0];
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

  // start kernel and allocate shared memory for "numberOfBins" floats
  cmtkDeviceHistogramEntropyKernel<<<dimGrid,dimBlock,numberOfBins*sizeof(float)>>>( result, dataPtr );
}

__global__
void 
cmtkDeviceHistogramPopulateWithMaskKernel( float* histPtr, const float *dataPtr, const int *maskPtr, const float rangeFrom, const float rangeTo, const int numberOfBins, const int iterations )
{
  int tx = threadIdx.x;

  // working histogram for this thread in shared memory
  float* working = &shared[(numberOfBins+1)*tx];

  // start by resetting all histogram bins to 0
  for ( int i = 0; i <= numberOfBins; ++i )
    working[i] = 0;

  // populate histogram bins
  const float binScale = (numberOfBins-1) / (rangeTo - rangeFrom);
  
  int offset = tx;
  for ( int i = 0; i < iterations; ++i, offset += blockDim.x )
    {
      float binIndex = fmaxf( 0, fminf( numberOfBins-1, (dataPtr[offset] - rangeFrom) * binScale ) );
      int index = truncf( (1+binIndex) * maskPtr[offset] );
      ++working[ index ];
    }

  // finally, add all thread working histograms to output histogram
  for ( int idx = 1+tx; idx <= numberOfBins; idx += blockDim.x )
    {
      float sum = 0;
      for ( int hx = 0; hx < blockDim.x; ++hx )
	{
	  sum += shared[ idx + hx*(1+numberOfBins) ];
	}
      histPtr[idx-1] += sum;
    }
}

__global__
void 
cmtkDeviceHistogramPopulateKernel( float* histPtr, const float *dataPtr, const float rangeFrom, const float rangeTo, const int numberOfBins, const int iterations )
{
  int tx = threadIdx.x;

  // working histogram for this thread in shared memory
  float* working = &shared[numberOfBins*tx];

  // start by resetting all histogram bins to 0
  for ( int i = 0; i < numberOfBins; ++i )
    working[i] = 0;

  // populate histogram bins
  const float binScale = (numberOfBins-1) / (rangeTo - rangeFrom);
  
  int offset = tx;
  for ( int i = 0; i < iterations; ++i, offset += blockDim.x )
    {
      int index = truncf( fmaxf( 0, fminf( numberOfBins-1, (dataPtr[offset] - rangeFrom) * binScale ) ) );
      ++working[ index ];
    }

  // finally, add all thread working histograms to output histogram
  for ( int idx = tx; idx < numberOfBins; idx += blockDim.x )
    {
      float sum = 0;
      for ( int hx = 0; hx < blockDim.x; ++hx )
	{
	  sum += shared[ idx + hx*numberOfBins ];
	}
      histPtr[idx] += sum;
    }
}

void
cmtkDeviceHistogramPopulate( float* histPtr, const float* dataPtr, const float rangeFrom, const float rangeTo, int numberOfBins, int numberOfSamples )
{
  // how many local copies of the histogram can we fit in shared memory?
  int device;
  cudaDeviceProp dprop;
  if ( (cudaGetDevice( &device ) != cudaSuccess) || (cudaGetDeviceProperties( &dprop, device ) != cudaSuccess ) )
    {
      fputs( "ERROR: could not get device properties.\n", stderr );
      exit( 1 );
    }
  
  int nThreads = dprop.sharedMemPerBlock / (sizeof(float) * numberOfBins);
  if ( nThreads > 512 )
    nThreads = 512;

  dim3 dimBlock( nThreads, 1 );
  dim3 dimGrid( 1, 1 );

  cmtkDeviceHistogramPopulateKernel<<<dimGrid,dimBlock,nThreads*numberOfBins*sizeof(float)>>>( histPtr, dataPtr, rangeFrom, rangeTo, numberOfBins, numberOfSamples / nThreads );

  const cudaError_t kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error code %d\n", static_cast<int>( kernelError ) );
      exit( 1 );      
    }

  const int residualSamples = numberOfSamples - nThreads * (numberOfSamples / nThreads);
  if ( residualSamples )
    {
      dim3 dimBlock( residualSamples, 1 );
      dim3 dimGrid( 1, 1 );
      
      cmtkDeviceHistogramPopulateKernel<<<dimGrid,dimBlock,residualSamples*numberOfBins*sizeof(float)>>>( histPtr, dataPtr + numberOfSamples - residualSamples, rangeFrom, rangeTo, numberOfBins, 1 );

      const cudaError_t kernelError = cudaGetLastError();
      if ( kernelError != cudaSuccess )
	{
	  fprintf( stderr, "ERROR: CUDA kernel failed with error code %d\n", static_cast<int>( kernelError ) );
	  exit( 1 );      
	}
    }
}

void
cmtkDeviceHistogramPopulate( float* histPtr, const float* dataPtr, const int* maskPtr, const float rangeFrom, const float rangeTo, int numberOfBins, int numberOfSamples )
{
  // how many local copies of the histogram can we fit in shared memory?
  int device;
  cudaDeviceProp dprop;
  if ( (cudaGetDevice( &device ) != cudaSuccess) || (cudaGetDeviceProperties( &dprop, device ) != cudaSuccess ) )
    {
      fputs( "ERROR: could not get device properties.\n", stderr );
      exit( 1 );
    }
  
  int nThreads = dprop.sharedMemPerBlock / (sizeof(float) * (1+numberOfBins));
  if ( nThreads > dprop.maxThreadsPerBlock )
    nThreads = dprop.maxThreadsPerBlock;

  dim3 dimBlock( nThreads, 1 );
  dim3 dimGrid( 1, 1 );
  
  cmtkDeviceHistogramPopulateWithMaskKernel<<<dimGrid,dimBlock,nThreads*(numberOfBins+1)*sizeof(float)>>>( histPtr, dataPtr, maskPtr, rangeFrom, rangeTo, numberOfBins, numberOfSamples / nThreads );
  
  const cudaError_t kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error code %d\n", static_cast<int>( kernelError ) );
      exit( 1 );      
    }
  
  const int residualSamples = numberOfSamples - nThreads * (numberOfSamples / nThreads);
  if ( residualSamples )
    {
      dim3 dimBlock( residualSamples, 1 );
      dim3 dimGrid( 1, 1 );
      
      cmtkDeviceHistogramPopulateWithMaskKernel<<<dimGrid,dimBlock,residualSamples*(numberOfBins+1)*sizeof(float)>>>( histPtr, dataPtr + numberOfSamples - residualSamples, maskPtr + numberOfSamples - residualSamples, rangeFrom, rangeTo, numberOfBins, 1 );

      const cudaError_t kernelError = cudaGetLastError();
      if ( kernelError != cudaSuccess )
	{
	  fprintf( stderr, "ERROR: CUDA kernel failed with error code %d\n", static_cast<int>( kernelError ) );
	  exit( 1 );      
	}
    }
}
