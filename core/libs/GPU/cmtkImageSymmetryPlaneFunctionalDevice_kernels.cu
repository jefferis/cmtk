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

#include "System/cmtkMemory.h"
#include "GPU/cmtkDeviceMemory.h"

#include <cuda_runtime_api.h>

#include <cstdio>

/// Texture reference to volume data.
texture<float, 3, cudaReadModeElementType> texRef;
texture<float, 3, cudaReadModeElementType> texRefX;

__global__
void cmtkImageSymmetryPlaneFunctionalDeviceEvaluateKernel( float* squares, const float matrix[4][4], const int dims0, const int dims1, const int dims2 )
{
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;

  const int z = blockIdx.x;

  const float mXz = z * matrix[2][0] + matrix[3][0];
  const float mYz = z * matrix[2][1] + matrix[3][1];
  const float mZz = z * matrix[2][2] + matrix[3][2];

  float sq = 0;
  for ( int y = ty; y < dims1; y += blockDim.y )
    {
      const float mXyz = y * matrix[1][0] + mXz;
      const float mYyz = y * matrix[1][1] + mYz;
      const float mZyz = y * matrix[1][2] + mZz;

      for ( int x = tx; x < dims0; x += blockDim.x )
	{
	  const float mX = x * matrix[0][0] + mXyz;
	  const float mY = x * matrix[0][1] + mYyz;
	  const float mZ = x * matrix[0][2] + mZyz;
	  
	  const float data = tex3D( texRef, x, y, z );
	  const float dataX = tex3D( texRefX, mX, mY, mZ );
	  
	  const float diff = data-dataX;
	  sq += diff*diff;
	}
    }

  const int idx = tx + blockDim.x * ( ty + blockDim.y * z );
  squares[idx] = sq;
}

__global__
void cmtkImageSymmetryPlaneFunctionalDeviceConsolidateKernel( float* squares, const int n )
{
  const int tx = threadIdx.x;

  for ( int i = tx + blockDim.x; i < n; i += blockDim.x )
    {
      squares[tx] += squares[i];
    }

  if ( tx == 0 )
    {
      for ( int i = 1; i < blockDim.x; ++i )
	squares[0] += squares[i];
    }
}

float
cmtkImageSymmetryPlaneFunctionalDeviceEvaluate( const int* dims3, void* array, const float matrix[4][4] )
{
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc( 32, 0, 0, 0, cudaChannelFormatKindFloat );
  
  // Set texture parameters for moving image interpolated access
  texRef.addressMode[0] = cudaAddressModeWrap;
  texRef.addressMode[1] = cudaAddressModeWrap;
  texRef.addressMode[2] = cudaAddressModeWrap;
  texRef.filterMode = cudaFilterModeLinear; 
  texRef.normalized = true; 

  // Bind the array to the texture reference 
  cudaError_t cudaError = cudaBindTextureToArray( texRef, (struct cudaArray*) array, channelDesc );
  if ( cudaError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: cudaBindTextureToArray failed with error '%s'\n", cudaGetErrorString( cudaError ) );
      exit( 1 );      
    }

  // Set texture parameters for fixed image indexed access
  texRefX.addressMode[0] = cudaAddressModeClamp;
  texRefX.addressMode[1] = cudaAddressModeClamp;
  texRefX.addressMode[2] = cudaAddressModeClamp;
  texRefX.filterMode = cudaFilterModePoint; 
  texRefX.normalized = false; 

  cudaError = cudaBindTextureToArray( texRefX, (struct cudaArray*) array, channelDesc );
  if ( cudaError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: cudaBindTextureToArray failed with error '%s'\n", cudaGetErrorString( cudaError ) );
      exit( 1 );      
    }

  // alocate memory for partial sums of squares
  cmtk::DeviceMemory<float>::SmartPtr partialSums = cmtk::DeviceMemory<float>::Create( 16*16*dims3[2] );

  dim3 dimBlock( 16, 16, 1 );
  dim3 dimGrid( dims3[2], 1 );
  
  cmtkImageSymmetryPlaneFunctionalDeviceEvaluateKernel<<<dimGrid,dimBlock>>>( partialSums->Ptr(), matrix, dims3[0], dims3[1], dims3[2] );

  cudaError = cudaGetLastError();
  if ( cudaError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error '%s'\n", cudaGetErrorString( cudaError ) );
      exit( 1 );      
    }

  const int nPixels = dims3[0]*dims3[1]*dims3[2];
  cmtkImageSymmetryPlaneFunctionalDeviceConsolidateKernel<<<dimGrid,dimBlock>>>( partialSums->Ptr(), nPixels );

  cudaError = cudaGetLastError();
  if ( cudaError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: CUDA kernel failed with error '%s'\n", cudaGetErrorString( cudaError ) );
      exit( 1 );      
    }

  cudaUnbindTexture( texRef );
  cudaUnbindTexture( texRefX );

  float result;
  partialSums->CopyToHost( &result, 1 );

  return result;
}
