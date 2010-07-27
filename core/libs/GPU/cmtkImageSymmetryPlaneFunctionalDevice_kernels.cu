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

#include "cmtkImageSymmetryPlaneFunctionalDevice_kernels.h"

#include <cuda_runtime_api.h>

#include <cstdio>

/// Texture reference to volume data.
texture<float, 3, cudaReadModeElementType> texRef;
texture<float, 3, cudaReadModeElementType> texRefX;

__global__
void cmtkImageSymmetryPlaneFunctionalDeviceEvaluateKernel( const float matrix[4][4], const float delta[3], const int dims0, const int dims1, const int dims2 )
{
  const int tx = threadIdx.x;

  __shared__ float sq[32];
  sq[tx] = 0;

  const int y = threadIdx.y;
  const int z = blockIdx.x;

  const float Y = y * delta[1];
  const float Z = z * delta[2];

  const float mXo = Y * matrix[1][0] + Z * matrix[2][0] + matrix[3][0];
  const float mYo = Y * matrix[1][1] + Z * matrix[2][1] + matrix[3][1];
  const float mZo = Y * matrix[1][2] + Z * matrix[2][2] + matrix[3][2];

  for ( int x = tx; x < dims0; x += blockDim.x )
    {
      const float X = x * delta[0];

      const float mX = X * matrix[0][0] + mXo;
      const float mY = X * matrix[0][1] + mYo;
      const float mZ = X * matrix[0][2] + mZo;

      const float data = tex3D( texRef, x, y, z );
      const float dataX = tex3D( texRefX, mX, mY, mZ );
      
      const float diff = data-dataX;
      sq[tx] += diff*diff;
    }

  // compute sum via butterfly
  for ( int bit = 1; bit < blockDim.x; bit <<= 1 )
    {
      const float sum = sq[tx] + sq[tx^bit];
      __syncthreads();
      sq[tx] = sum;
      __syncthreads();
    }
}

float
cmtkImageSymmetryPlaneFunctionalDeviceEvaluate( const int* dims3, void* array, const float matrix[4][4], const float delta[3] )
{
  // Set texture parameters
  texRef.addressMode[0] = cudaAddressModeWrap;
  texRef.addressMode[1] = cudaAddressModeWrap;
  texRef.addressMode[2] = cudaAddressModeWrap;
  texRef.filterMode = cudaFilterModeLinear; 
  texRef.normalized = true; 

  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc( 32, 0, 0, 0, cudaChannelFormatKindFloat );
  
  // Bind the array to the texture reference 
  cudaBindTextureToArray( texRef, (struct cudaArray*) array, channelDesc );

  // Set texture parameters
  texRefX.addressMode[0] = cudaAddressModeClamp;
  texRefX.addressMode[1] = cudaAddressModeClamp;
  texRefX.addressMode[2] = cudaAddressModeClamp;
  texRefX.filterMode = cudaFilterModePoint; 
  texRefX.normalized = false; 

  cudaBindTextureToArray( texRefX, (struct cudaArray*) array, channelDesc );

  dim3 dimBlock( 32, dims3[1], 1 );
  dim3 dimGrid( dims3[2], 1 );
  
  cmtkImageSymmetryPlaneFunctionalDeviceEvaluateKernel<<<dimGrid,dimBlock>>>( matrix, delta, dims3[0], dims3[1], dims3[2] );

  const cudaError_t kernelError = cudaGetLastError();
  if ( kernelError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error %s\n",cudaGetErrorString( kernelError ) );
      exit( 1 );      
    }

  cudaUnbindTexture( texRef );

  return 0;
}
