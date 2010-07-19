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
cmtkDeviceHistogramConsolidate( float* histPtr, float* localPtr, const int numberOfBins, const int numberOfThreads )
{
  int tx = threadIdx.x;

  // finally, add all thread working histograms to output histogram
  for ( int idx = 1+tx; idx <= numberOfBins; idx += blockDim.x )
    {
      float sum = 0;
      for ( int hx = 0; hx < numberOfThreads; ++hx )
	{
	  sum += localPtr[ idx + hx*(1+numberOfBins) ];
	}
      histPtr[idx-1] += sum;
    }
}

__global__
void 
cmtkDeviceHistogramPopulateKernel( float* histPtr, float* localPtr, const float *dataPtr, const float rangeFrom, const float rangeTo, const int numberOfBins, const int numberOfSamples )
{
  int tx = threadIdx.x;
  int offs = tx + blockDim.x * blockIdx.x;
  int skip = blockDim.x * gridDim.x;

  // working histogram for this thread
  float* working = localPtr + (numberOfBins+1)*offs;

  // start by resetting all histogram bins to 0
  for ( int i = 0; i < numberOfBins; ++i )
    working[i] = 0;

  // populate histogram bins
  const float binScale = (numberOfBins-1) / (rangeTo - rangeFrom);
  
  for ( int offset = offs; offset < numberOfSamples; offset += skip )
    {
      int index = 1+truncf( fmaxf( 0, fminf( numberOfBins-1, (dataPtr[offset] - rangeFrom) * binScale ) ) ); // 1+... for consistency with masked computation; bin0 is ignored in final analysis.
      ++working[ index ];
    }
}

__global__
void 
cmtkDeviceHistogramPopulateLogKernel( float* histPtr, float* localPtr, const float *dataPtr, const float rangeFrom, const float rangeTo, const int numberOfBins, const int numberOfSamples )
{
  int tx = threadIdx.x;
  int offs = tx + blockDim.x * blockIdx.x;
  int skip = blockDim.x * gridDim.x;

  // working histogram for this thread
  float* working = localPtr + (numberOfBins+1)*offs;

  // start by resetting all histogram bins to 0
  for ( int i = 0; i < numberOfBins; ++i )
    working[i] = 0;

  // populate histogram bins
  const float binScale = (numberOfBins-1) / (rangeTo - rangeFrom);
  const float logNumBins = log( static_cast<float>( numberOfBins ) );
  
  for ( int offset = offs; offset < numberOfSamples; offset += skip )
    {
      int index = 1+truncf( (numberOfBins-1) * fmaxf( 0, fminf( 1, log((1 + dataPtr[offset]-rangeFrom)*binScale)/logNumBins ) ) );
      ++working[ index ];
    }
}

void
cmtkDeviceHistogramPopulate( float* histPtr, const float* dataPtr, const float rangeFrom, const float rangeTo, const bool logScale, const int numberOfBins, const int numberOfSamples )
{
  dim3 dimBlock( 256, 1 );
  dim3 dimGrid( 16, 1 );

  const size_t nThreads = 16*256;
  const size_t lBytes = sizeof(float) * (numberOfBins+1) * nThreads;

  float* localPtr;
  if ( cudaMalloc( &localPtr, lBytes ) != cudaSuccess )
    {
      fprintf( stderr, "ERROR: cudaMalloc() failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );
    }

  if ( cudaMemset( localPtr, 0, lBytes ) != cudaSuccess )
    {
      fprintf( stderr, "ERROR: cudaMemset() failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );
    }

  if ( logScale )
    cmtkDeviceHistogramPopulateLogKernel<<<dimGrid,dimBlock>>>( histPtr, localPtr, dataPtr, rangeFrom, rangeTo, numberOfBins, numberOfSamples );
  else
    cmtkDeviceHistogramPopulateKernel<<<dimGrid,dimBlock>>>( histPtr, localPtr, dataPtr, rangeFrom, rangeTo, numberOfBins, numberOfSamples );

  cudaError_t kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error %s\n", cudaGetErrorString( kernelError ) );
      exit( 1 );      
    }

  dim3 dimBlock2( 256, 1 );
  dim3 dimGrid2( 1, 1 );

  cmtkDeviceHistogramConsolidate<<<dimGrid2,dimBlock2>>>( histPtr, localPtr, numberOfBins, nThreads );
  
  kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error %s\n", cudaGetErrorString( kernelError ) );
      exit( 1 );      
    }

  cudaFree( localPtr );
}


__global__
void 
cmtkDeviceHistogramPopulateWithMaskKernel( float* histPtr, float* localPtr, const float *dataPtr, const int *maskPtr, const float rangeFrom, const float rangeTo, const int numberOfBins, const int numberOfSamples )
{
  int tx = threadIdx.x;
  int offs = tx + blockDim.x * blockIdx.x;
  int skip = blockDim.x * gridDim.x;

  // working histogram for this thread in shared memory
  float* working = localPtr + (numberOfBins+1)*offs;

  // populate histogram bins
  const float binScale = (numberOfBins-1) / (rangeTo - rangeFrom);
  
  for ( int offset = offs; offset < numberOfSamples; offset += skip )
    {
      const float d = (dataPtr[offset] - rangeFrom) * binScale;
      const float m = maskPtr[offset];
      
      float binIndex = fmaxf( 0, fminf( numberOfBins-1, d ) );
      const int index = truncf( (1+binIndex) * m );
      
      ++working[ index ];
    }
}

__global__
void 
cmtkDeviceHistogramPopulateLogWithMaskKernel( float* histPtr, float* localPtr, const float *dataPtr, const int *maskPtr, const float rangeFrom, const float rangeTo, const int numberOfBins, const int numberOfSamples )
{
  int tx = threadIdx.x;
  int offs = tx + blockDim.x * blockIdx.x;
  int skip = blockDim.x * gridDim.x;

  // working histogram for this thread in shared memory
  float* working = localPtr + (numberOfBins+1)*offs;

  // populate histogram bins
  const float binScale = (numberOfBins-1) / (rangeTo - rangeFrom);
  const float logNumBins = log( static_cast<float>( numberOfBins ) );
  
  for ( int offset = offs; offset < numberOfSamples; offset += skip )
    {
      const float d = log((1 + dataPtr[offset]-rangeFrom)*binScale);
      const float m = maskPtr[offset];
      
      float binIndex = (numberOfBins-1) * fmaxf( 0, fminf( 1, d / logNumBins ) );
      const int index = truncf( (1+binIndex) * m );
      ++working[ index ];
    }
}

void
cmtkDeviceHistogramPopulate( float* histPtr, const float* dataPtr, const int* maskPtr, const float rangeFrom, const float rangeTo, const bool logScale, const int numberOfBins, const int numberOfSamples )
{
  dim3 dimBlock( 256, 1 );
  dim3 dimGrid( 16, 1 );

  const size_t nThreads = 16*256;
  const size_t lBytes = sizeof(float) * (numberOfBins+1) * nThreads;

  float* localPtr;
  if ( cudaMalloc( &localPtr, lBytes ) != cudaSuccess )
    {
      fprintf( stderr, "ERROR: cudaMalloc() failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );
    }

  if ( cudaMemset( localPtr, 0, lBytes ) != cudaSuccess )
    {
      fprintf( stderr, "ERROR: cudaMemset() failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );
    }

  if ( logScale )
    cmtkDeviceHistogramPopulateLogWithMaskKernel<<<dimGrid,dimBlock>>>( histPtr, localPtr, dataPtr, maskPtr, rangeFrom, rangeTo, numberOfBins, numberOfSamples );
  else
    cmtkDeviceHistogramPopulateWithMaskKernel<<<dimGrid,dimBlock>>>( histPtr, localPtr, dataPtr, maskPtr, rangeFrom, rangeTo, numberOfBins, numberOfSamples );

  cudaError_t kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error %s\n", cudaGetErrorString( kernelError ) );
      exit( 1 );      
    }

  dim3 dimBlock2( 256, 1 );
  dim3 dimGrid2( 1, 1 );

  cmtkDeviceHistogramConsolidate<<<dimGrid2,dimBlock2>>>( histPtr, localPtr, numberOfBins, nThreads );
  
  kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error %s\n", cudaGetErrorString( kernelError ) );
      exit( 1 );      
    }

  cudaFree( localPtr );
}
