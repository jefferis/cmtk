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

#include "GPU/cmtkCUDA.h"
#include "GPU/cmtkDeviceMemory.h"
#include "GPU/cmtkSumReduction_kernel.h"

#include <cuda_runtime_api.h>

/// Texture reference to volume data.
texture<float, 3, cudaReadModeElementType> texRefMov;
texture<float, 3, cudaReadModeElementType> texRefFix;

__constant__ float deviceMatrix[4][4];

__global__
void
cmtkImageSymmetryPlaneFunctionalDeviceEvaluateMSDKernel( float* squares, const int dims0, const int dims1, const int dims2 )
{
  const int tx = threadIdx.x;
  const int ty = threadIdx.y;

  float sq = 0;

  for ( int z = blockIdx.x; z < dims2; z += blockDim.x )
    {
      const float mXz = z * deviceMatrix[2][0] + deviceMatrix[3][0];
      const float mYz = z * deviceMatrix[2][1] + deviceMatrix[3][1];
      const float mZz = z * deviceMatrix[2][2] + deviceMatrix[3][2];
      
      for ( int y = ty; y < dims1; y += blockDim.y )
	{
	  const float mXyz = y * deviceMatrix[1][0] + mXz;
	  const float mYyz = y * deviceMatrix[1][1] + mYz;
	  const float mZyz = y * deviceMatrix[1][2] + mZz;
	  
	  for ( int x = tx; x < dims0; x += blockDim.x )
	    {
	      const float dataFix = tex3D( texRefFix, x, y, z );

	      const float mX = x * deviceMatrix[0][0] + mXyz;
	      const float mY = x * deviceMatrix[0][1] + mYyz;
	      const float mZ = x * deviceMatrix[0][2] + mZyz;
	      
	      const float dataMov = tex3D( texRefMov, mX, mY, mZ );
	      
	      const float diff = dataMov-dataFix;
	      sq += diff*diff;
	      sq += deviceMatrix[0][0];
	    }
	}
    }

  const int idx = tx + blockDim.x * ( ty + blockDim.y * blockIdx.x );
  squares[idx] = sq;
}

float
cmtk::ImageSymmetryPlaneFunctionalDeviceEvaluateMSD( const int* dims3, void* array, const float matrix[4][4] )
{
  cmtkCheckCallCUDA( cudaMemcpyToSymbol( deviceMatrix, matrix, 16 * sizeof( float ), 0, cudaMemcpyHostToDevice ) );

  // Set texture parameters for fixed image indexed access
  texRefFix.addressMode[0] = cudaAddressModeClamp;
  texRefFix.addressMode[1] = cudaAddressModeClamp;
  texRefFix.addressMode[2] = cudaAddressModeClamp;
  texRefFix.filterMode = cudaFilterModePoint; 
  texRefFix.normalized = false; 

  cmtkCheckCallCUDA( cudaBindTextureToArray( texRefFix, (struct cudaArray*) array, cudaCreateChannelDesc<float>() ) );

  // Set texture parameters for moving image interpolated access
  texRefMov.addressMode[0] = cudaAddressModeWrap;
  texRefMov.addressMode[1] = cudaAddressModeWrap;
  texRefMov.addressMode[2] = cudaAddressModeWrap;
  texRefMov.filterMode = cudaFilterModeLinear; 
  texRefMov.normalized = true; 

  // Bind the array to the texture reference 
  cmtkCheckCallCUDA( cudaBindTextureToArray( texRefMov, (struct cudaArray*) array, cudaCreateChannelDesc<float>() ) );

  // alocate memory for partial sums of squares
  dim3 dimBlock( 16, 16 );

  const int nPartials = dimBlock.x * dimBlock.y * dimBlock.z * dims3[2];
  cmtk::DeviceMemory<float> partialSums( nPartials );
  
  cmtkImageSymmetryPlaneFunctionalDeviceEvaluateMSDKernel<<<dims3[2],dimBlock>>>( partialSums.Ptr(), dims3[0], dims3[1], dims3[2] );
  cmtkCheckLastErrorCUDA;
  
  cmtkCheckCallCUDA( cudaUnbindTexture( texRefMov ) );
  cmtkCheckCallCUDA( cudaUnbindTexture( texRefFix ) );

  return cmtk::SumReduction( partialSums.Ptr(), nPartials ) / dims3[0]*dims3[1]*dims3[2];
}
