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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_kernels.h"

#include <cstdio>
#include <algorithm>

__constant__ float deviceWeights[19];
__constant__ float deviceCorrections[19];

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalUpdateOutputImageKernel
( float* output, float* input, int degree, int multiply, int nPixels, int dims0, int dims1, int dims2 )
{  
  const int offset = blockIdx.x * blockDim.x + threadIdx.x;

  if ( offset < nPixels )
    {
      const int z = offset / (dims1*dims2);
      const int y = (offset - z*(dims1*dims2)) / dims1;
      const int x = offset % dims1;

      const float X = 2.0f * (x-dims0/2) / dims0;
      const float Y = 2.0f * (y-dims1/2) / dims1;
      const float Z = 2.0f * (z-dims2/2) / dims2;
      
      const float in = input[offset];
      
      float bias =
	deviceWeights[0] * (X - deviceCorrections[0]) + 
	deviceWeights[1] * (Y - deviceCorrections[1]) +
	deviceWeights[2] * (Z - deviceCorrections[2]);
      
      if ( degree > 1 )
	{
	  bias +=
	    deviceWeights[3] * (X * X - deviceCorrections[3]) +
	    deviceWeights[4] * (X * Y - deviceCorrections[4]) +
	    deviceWeights[5] * (X * Z - deviceCorrections[5]) +
	    deviceWeights[6] * (Y * Y - deviceCorrections[6]) +
	    deviceWeights[7] * (Y * Z - deviceCorrections[7]) +
	    deviceWeights[8] * (Z * Z - deviceCorrections[8]);
	}
      
      if ( degree > 2 )
	{
	  bias +=
	    deviceWeights[ 9] * (X * X * X - deviceCorrections[ 9]) +
	    deviceWeights[10] * (X * X * Y - deviceCorrections[10]) +
	    deviceWeights[11] * (X * X * Z - deviceCorrections[11]) +
	    deviceWeights[12] * (X * Y * Y - deviceCorrections[12]) +
	    deviceWeights[13] * (X * Y * Z - deviceCorrections[13]) +
	    deviceWeights[14] * (X * Z * Z - deviceCorrections[14]) +
	    deviceWeights[15] * (Y * Y * Y - deviceCorrections[15]) +
	    deviceWeights[16] * (Y * Y * Z - deviceCorrections[16]) +
	    deviceWeights[17] * (Y * Z * Z - deviceCorrections[17]) +
	    deviceWeights[18] * (Z * Z * Z - deviceCorrections[18]);
	}

      if ( multiply )
	{
	  output[offset] = in * (bias+1);
	}
      else
	{
	  output[offset] = in + bias;
	}
    }
}

void
cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage
( float* output, float* input, const int dims0, const int dims1, const int dims2, const int degree, const int multiply, const int nargs, const float* weights, const float* corrections )
{ 
  const int nPixels = dims0 * dims1 * dims2;

  // how many local copies of the histogram can we fit in shared memory?
  int device;
  cudaDeviceProp dprop;
  if ( (cudaGetDevice( &device ) != cudaSuccess) || (cudaGetDeviceProperties( &dprop, device ) != cudaSuccess ) )
    {
      fprintf( stderr, "ERROR: cudaGetDevice() failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );
    }
  
  const int nThreads = std::min<int>( nPixels, dprop.maxThreadsPerBlock );
  
  dim3 dimBlock( nThreads, 1, 1 );
  dim3 dimGrid( 1+(nPixels-1)/nThreads, 1 );
  
  if ( (cudaMemcpy( deviceWeights, weights, nargs * sizeof( *weights ), cudaMemcpyHostToDevice ) != cudaSuccess) ||
       (cudaMemcpy( deviceCorrections, corrections, nargs * sizeof( *corrections ), cudaMemcpyHostToDevice ) != cudaSuccess) )
    {
      fprintf( stderr, "ERROR: cudaMemcpy() to constant memory failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );      
    }
    
  cmtkEntropyMinimizationIntensityCorrectionFunctionalUpdateOutputImageKernel<<<dimGrid,dimBlock>>>( output, input, degree, multiply, nPixels, dims0, dims1, dims2 );

  const cudaError_t kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error %s\n",cudaGetErrorString( kernelError ) );
      exit( 1 );      
    }
}
