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

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageAddKernel
( float* output, float* input, float* add )
{
  int x = blockIdx.x * 512 + threadIdx.x;
  output[x] = input[x] + add[x];
}

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageMulKernel
( float* output, float* input, float* mul )
{
  int x = blockIdx.x * 512 + threadIdx.x;
  output[x] = input[x] * mul[x];
}

void
cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( float* input, float* output, float* biasAdd, float* biasMul, int numberOfPixels )
{
  dim3 dimBlock( 512, 1 );
  dim3 dimGrid( 1+((numberOfPixels-1)/512), 1 );

  if ( biasMul )
    cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageMulKernel<<<dimGrid,dimBlock>>>( output, input, biasMul );

  if ( biasAdd )
    cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageAddKernel<<<dimGrid,dimBlock>>>( output, input, biasAdd );
}

__constant__ float weights[19];
__constant__ float correction[19];

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAComputeMonomials
( float* output, int slice, int dims0, int dims1, int dims2 )
{
  int x = blockIdx.x * 16 + threadIdx.x;
  int y = blockIdx.y * 16 + threadIdx.y;
  int z = threadIdx.z + slice;

  float X = 2.0 * (x-dims0/2) / dims0;
  float Y = 2.0 * (y-dims1/2) / dims1;
  float Z = 2.0 * (z-dims2/2) / dims2;

  int offset = x + dims0 * (y + dims1 * (z+slice) );
  float bias =
    weights[0] * (X - correction[0]) + 
    weights[1] * (Y - correction[1]) +
    weights[2] * (Z - correction[2]);

  bias +=
    weights[3] * (X * X - correction[3])+
    weights[4] * (X * Y - correction[4])+
    weights[5] * (X * Z - correction[5]) +
    weights[6] * (Y * Y - correction[6]) +
    weights[7] * (Y * Z - correction[7]) +
    weights[8] * (Z * Z - correction[8]);
  
  bias +=
    weights[ 9] * (X * X * X - correction[ 9]) +
    weights[10] * (X * X * Y - correction[10]) +
    weights[11] * (X * X * Z - correction[11]) +
    weights[12] * (X * Y * Y - correction[12]) +
    weights[13] * (X * Y * Z - correction[13]) +
    weights[14] * (X * Z * Z - correction[14]) +
    weights[15] * (Y * Y * Y - correction[15]) +
    weights[16] * (Y * Y * Z - correction[16]) +
    weights[17] * (Y * Z * Z - correction[17]) +
    weights[18] * (Z * Z * Z - correction[18]);
}
