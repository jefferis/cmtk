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

__constant__ float deviceWeights[19];
__constant__ float deviceCorrections[19];

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalUpdateOutputImageKernel
( float* output, float* input, int degree, int multiply, int slice, int dims0, int dims1, int dims2 )
{
  int x = blockIdx.x * 16 + threadIdx.x;
  int y = blockIdx.y * 16 + threadIdx.y;
  int z = threadIdx.z + slice;

  if ( (x >= dims0) || (y >= dims1) || (z >= dims2) )
    return;

  float X = 2.0 * (x-dims0/2) / dims0;
  float Y = 2.0 * (y-dims1/2) / dims1;
  float Z = 2.0 * (z-dims2/2) / dims2;

  int offset = x + dims0 * (y + dims1 * (z+slice) );
  const float in = input[offset];
  
  float bias =
    deviceWeights[0] * (X - deviceCorrections[0]) + 
    deviceWeights[1] * (Y - deviceCorrections[1]) +
    deviceWeights[2] * (Z - deviceCorrections[2]);

  if ( degree > 1 )
    {
      bias +=
	deviceWeights[3] * (X * X - deviceCorrections[3])+
	deviceWeights[4] * (X * Y - deviceCorrections[4])+
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
      output[offset] = in * bias;
    }
  else
    {
      output[offset] = in + bias;
    }
}

void
cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage
( float* output, float* input, const int dims0, const int dims1, const int dims2, const int degree, const int multiply, const int nargs, const float* weights, const float* corrections )
{ 
  const int planesPerSlice = 8; // create 16*16*8 = 512 threads per block
  dim3 dimBlock( 16, 16, planesPerSlice );

  const int nSlices = 1+((dims2-1)/planesPerSlice);
  dim3 dimGrid( 1+((dims0-1)/16), 1+((dims1-1)/16), nSlices );

  cudaMemcpy( deviceWeights, weights, nargs * sizeof( *weights ), cudaMemcpyHostToDevice );
  cudaMemcpy( deviceCorrections, corrections, nargs * sizeof( *corrections ), cudaMemcpyHostToDevice );
  
  for ( int slice = 0; slice < dims2; slice += planesPerSlice )
    {
      cmtkEntropyMinimizationIntensityCorrectionFunctionalUpdateOutputImageKernel<<<dimGrid,dimBlock>>>( output, input, degree, multiply, slice, dims0, dims1, dims2 );
    }
}
