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

#include "GPU/cmtkCUDA.h"
#include "GPU/cmtkDeviceMemory.h"

#include <cuda_runtime_api.h>

/// Texture reference to volume data.
texture<float, 3, cudaReadModeElementType> texRef;

__constant__ float deviceKernel[128];

__global__
void
cmtkDeviceImageConvolutionKernelX( float* dest, int dims0, int dims1, int dims2, int kernelLength, int kernelCenter )
{
  const int offs = threadIdx.x + blockDim.x * blockIdx.x;

  const int x = offs % dims0;
  const int y = (offs / dims0) % dims1;
  const int z = offs / (dims0 * dims1);

  float sum = 0, total = 0;
  for ( int i = 0; i < kernelLength; ++i )
    {
      const int xx = x + i - kernelCenter;
      
      const float w = ( (xx>=0) && (xx<dims0) ) ? deviceKernel[i] : 0;
      
      sum += w * tex3D( texRef, xx, y, z );
      total += w;
    }
  
  dest[offs] = sum / total;
}

__global__
void
cmtkDeviceImageConvolutionKernelY( float* dest, int dims0, int dims1, int dims2, int kernelLength, int kernelCenter )
{
  const int offs = threadIdx.x + blockDim.x * blockIdx.x;

  const int x = offs % dims0;
  const int y = (offs / dims0) % dims1;
  const int z = offs / (dims0 * dims1);

  float sum = 0, total = 0;
  for ( int i = 0; i < kernelLength; ++i )
    {
      const int yy = y + i - kernelCenter;
      
      const float w = ( (yy>=0) && (yy<dims1) ) ? deviceKernel[i] : 0;
      
      sum += w * tex3D( texRef, x, yy, z );
      total += w;
    }
  
  dest[offs] = sum / total;
}

__global__
void
cmtkDeviceImageConvolutionKernelZ( float* dest, int dims0, int dims1, int dims2, int kernelLength, int kernelCenter )
{
  const int offs = threadIdx.x + blockDim.x * blockIdx.x;

  const int x = offs % dims0;
  const int y = (offs / dims0) % dims1;
  const int z = offs / (dims0 * dims1);

  float sum = 0, total = 0;
  for ( int i = 0; i < kernelLength; ++i )
    {
      const int zz = z + i - kernelCenter;
      
      const float w = ( (zz>=0) && (zz<dims2) ) ? deviceKernel[i] : 0;
      
      sum += w * tex3D( texRef, x, y, zz );
      total += w;
    }
  
  dest[offs] = sum / total;
}

void
cmtk::DeviceImageConvolution( float* dest, const int* dims3, void* array, const int kernelLengthX, const float* kernelX, const int kernelLengthY, const float* kernelY, const int kernelLengthZ, const float* kernelZ )
{
  // Set texture parameters for fixed image indexed access
  texRef.addressMode[0] = cudaAddressModeClamp;
  texRef.addressMode[1] = cudaAddressModeClamp;
  texRef.addressMode[2] = cudaAddressModeClamp;
  texRef.filterMode = cudaFilterModePoint; 
  texRef.normalized = false; 

  cmtkCheckCallCUDA( cudaBindTextureToArray( texRef, (struct cudaArray*) array, cudaCreateChannelDesc<float>() ) );

  const int nPixels = dims3[0] * dims3[1] * dims3[2];

  cmtkCheckCallCUDA( cudaMemcpyToSymbol( deviceKernel, kernelX, kernelLengthX * sizeof( float ), 0, cudaMemcpyHostToDevice ) );
  
  dim3 threads( 32 );
  dim3 blocks( nPixels/32+1 );

  cmtkDeviceImageConvolutionKernelX<<<threads,blocks>>>( dest, dims3[0], dims3[1], dims3[2], kernelLengthX, (kernelLengthX-1)>>1 );
  cmtkCheckLastErrorCUDA;

  cmtkCheckCallCUDA( cudaMemcpyToArray( (struct cudaArray*) array, 0, 0, dest, nPixels, cudaMemcpyDeviceToDevice ) );

  cmtkCheckCallCUDA( cudaMemcpyToSymbol( deviceKernel, kernelY, kernelLengthY * sizeof( float ), 0, cudaMemcpyHostToDevice ) );
  cmtkDeviceImageConvolutionKernelY<<<threads,blocks>>>( dest, dims3[0], dims3[1], dims3[2], kernelLengthY, (kernelLengthY-1)>>1 );
  cmtkCheckLastErrorCUDA;

  cmtkCheckCallCUDA( cudaMemcpyToArray( (struct cudaArray*) array, 0, 0, dest, nPixels, cudaMemcpyDeviceToDevice ) );

  cmtkCheckCallCUDA( cudaMemcpyToSymbol( deviceKernel, kernelZ, kernelLengthZ * sizeof( float ), 0, cudaMemcpyHostToDevice ) );  
  cmtkDeviceImageConvolutionKernelZ<<<threads,blocks>>>( dest, dims3[0], dims3[1], dims3[2], kernelLengthZ, (kernelLengthZ-1)>>1 );
  cmtkCheckLastErrorCUDA;
  
  cmtkCheckCallCUDA( cudaUnbindTexture( texRef ) );
}

void
cmtk::DeviceImageConvolutionInPlace( const int* dims3, void* array, const int kernelLengthX, const float* kernelX, const int kernelLengthY, const float* kernelY, const int kernelLengthZ, const float* kernelZ )
{
  const int nPixels = dims3[0] * dims3[1] * dims3[2];
  cmtk::DeviceMemory<float>::SmartPtr temporary = cmtk::DeviceMemory<float>::Create( nPixels );

  // call out-of-place-place convolution

  DeviceImageConvolution( temporary->Ptr(), dims3, array, kernelLengthX, kernelX, kernelLengthY, kernelY, kernelLengthZ, kernelZ );
  
  // copy back into original array
  cmtkCheckCallCUDA( cudaMemcpyToArray( (struct cudaArray*) array, 0, 0, temporary->Ptr(), nPixels, cudaMemcpyDeviceToDevice ) );
  
  cmtkCheckCallCUDA( cudaUnbindTexture( texRef ) );
}
