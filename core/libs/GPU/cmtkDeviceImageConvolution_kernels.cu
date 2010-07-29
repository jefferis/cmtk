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

#include "cmtkDeviceImageConvolution_kernels.h"

#include "System/cmtkMemory.h"
#include "GPU/cmtkDeviceMemory.h"

#include <cuda_runtime_api.h>

#include <cstdio>

/// Texture reference to volume data.
texture<float, 3, cudaReadModeElementType> texRef;

__constant__ float deviceKernel[128];

__global__
void
cmtkDeviceImageConvolutionKernelX( float* dest, int dims0, int dims1, int dims2, int kernelLength )
{
}

__global__
void
cmtkDeviceImageConvolutionKernelY( float* dest, int dims0, int dims1, int dims2, int kernelLength )
{
}

__global__
void
cmtkDeviceImageConvolutionKernelZ( float* dest, int dims0, int dims1, int dims2, int kernelLength )
{
}

void
cmtkDeviceImageConvolution( const int* dims3, void* array, const int kernelLengthX, const float* kernelX, const int kernelLengthY, const float* kernelY, const int kernelLengthZ, const float* kernelZ )
{
  // Set texture parameters for fixed image indexed access
  texRef.addressMode[0] = cudaAddressModeClamp;
  texRef.addressMode[1] = cudaAddressModeClamp;
  texRef.addressMode[2] = cudaAddressModeClamp;
  texRef.filterMode = cudaFilterModePoint; 
  texRef.normalized = false; 

  cudaError_t cudaError = cudaBindTextureToArray( texRef, (struct cudaArray*) array, cudaCreateChannelDesc<float>() );
  if ( cudaError != cudaSuccess )
    {
      fprintf( stderr, "ERROR: cudaBindTextureToArray failed with error '%s'\n", cudaGetErrorString( cudaError ) );
      exit( 1 );      
    }

  const int nPixels = dims3[0] * dims3[1] * dims3[2];
  cmtk::DeviceMemory<float>::SmartPtr temporary = cmtk::DeviceMemory<float>::Create( nPixels );

  if ( (cudaMemcpyToSymbol( deviceKernel, kernelX, kernelLengthX * sizeof( float ), 0, cudaMemcpyHostToDevice ) != cudaSuccess) )
    {
      fprintf( stderr, "ERROR: cudaMemcpyToSymbol() to constant memory failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );      
    }
  
  dim3 threads;
  dim3 blocks;

  cmtkDeviceImageConvolutionKernelX<<<threads,blocks>>>( temporary->Ptr(), dims3[0], dims3[1], dims3[2], kernelLengthX );

  cudaMemcpyToArray( (struct cudaArray*) array, 0, 0, temporary->Ptr(), nPixels, cudaMemcpyDeviceToDevice );

  if ( (cudaMemcpyToSymbol( deviceKernel, kernelY, kernelLengthY * sizeof( float ), 0, cudaMemcpyHostToDevice ) != cudaSuccess) )
    {
      fprintf( stderr, "ERROR: cudaMemcpyToSymbol() to constant memory failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );      
    }
  
  cmtkDeviceImageConvolutionKernelY<<<threads,blocks>>>( temporary->Ptr(), dims3[0], dims3[1], dims3[2], kernelLengthY );

  cudaMemcpyToArray( (struct cudaArray*) array, 0, 0, temporary->Ptr(), nPixels, cudaMemcpyDeviceToDevice );

  if ( (cudaMemcpyToSymbol( deviceKernel, kernelZ, kernelLengthZ * sizeof( float ), 0, cudaMemcpyHostToDevice ) != cudaSuccess) )
    {
      fprintf( stderr, "ERROR: cudaMemcpyToSymbol() to constant memory failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );      
    }
  
  cmtkDeviceImageConvolutionKernelZ<<<threads,blocks>>>( temporary->Ptr(), dims3[0], dims3[1], dims3[2], kernelLengthZ );

  cudaMemcpyToArray( (struct cudaArray*) array, 0, 0, temporary->Ptr(), nPixels, cudaMemcpyDeviceToDevice );
  
  cudaUnbindTexture( texRef );
}
