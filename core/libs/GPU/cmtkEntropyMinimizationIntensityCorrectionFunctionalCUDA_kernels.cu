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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDA_kernels.h"

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
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImage( float* input, float* output, float* biasAdd, float* biasMul, int numberOfPixels )
{
  dim3 dimBlock( 512, 1 );
  dim3 dimGrid( 1+((numberOfPixels-1)/512), 1 );

  if ( biasMul )
    cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageMulKernel<<<dimGrid,dimBlock>>>( output, input, biasMul );

  if ( biasAdd )
    cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageAddKernel<<<dimGrid,dimBlock>>>( output, input, biasAdd );
}

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAComputeMonomials1
( float* output, float* weights, int slice, int dims0, int dims1, int dims2 )
{
  int x = blockIdx.x * 16 + threadIdx.x;
  int y = blockIdx.y * 16 + threadIdx.y;
  int z = threadIdx.z + slice;

  int offset = x + dims0 * (y + dims1 * z );
  output[offset] +=
    weights[0] * 2.0 * (x-dims0/2) / dims0 + 
    weights[1] * 2.0 * (y-dims1/2) / dims1 +
    weights[2] * 2.0 * (z-dims2/2) / dims2;
}

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAComputeMonomials2
( float* output, float* weights, int slice, int dims0, int dims1, int dims2 )
{
  int x = blockIdx.x * 16 + threadIdx.x;
  int y = blockIdx.y * 16 + threadIdx.y;
  int z = threadIdx.z + slice;

  float X = 2.0 * (x-dims0/2) / dims0;
  float Y = 2.0 * (y-dims1/2) / dims1;
  float Z = 2.0 * (z-dims2/2) / dims2;

  int offset = x + dims0 * (y + dims1 * (z+slice) );
  output[offset] +=
    weights[0] * X * X +
    weights[1] * X * Y +
    weights[2] * X * Z +
    weights[3] * Y * Y +
    weights[4] * Y * Z +
    weights[5] * Z * Z;
}

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAComputeMonomials3
( float* output, float* weights, int slice, int dims0, int dims1, int dims2 )
{
  int x = blockIdx.x * 16 + threadIdx.x;
  int y = blockIdx.y * 16 + threadIdx.y;
  int z = threadIdx.z + slice;

  float X = 2.0 * (x-dims0/2) / dims0;
  float Y = 2.0 * (y-dims1/2) / dims1;
  float Z = 2.0 * (z-dims2/2) / dims2;

  int offset = x + dims0 * (y + dims1 * (z+slice) );
  output[offset] +=
    weights[ 0] * X * X * X +
    weights[ 1] * X * X * Y +
    weights[ 2] * X * X * Z +
    weights[ 3] * X * Y * Y +
    weights[ 4] * X * Y * Z +
    weights[ 5] * X * Z * Z +
    weights[ 6] * Y * Y * Y +
    weights[ 7] * Y * Y * Z +
    weights[ 8] * Y * Z * Z +
    weights[ 9] * Z * Z * Z;
}
