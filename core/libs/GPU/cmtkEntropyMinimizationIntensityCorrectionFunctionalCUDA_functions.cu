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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDA_functions.h"

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageAddKernel
( float* output, float* input, float* add, int numberOfPixels )
{
  int tx = threadIdx.x;

  output[tx] = input[tx] + add[tx];
}

__global__
void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageMulKernel
( float* output, float* input, float* mul, int numberOfPixels )
{
  int tx = threadIdx.x;

  output[tx] = input[tx] * mul[tx];
}

void
cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImage( float* input, float* output, float* biasAdd, float* biasMul, int numberOfPixels )
{
  dim3 dimBlock( 512, 1 );
  dim3 dimGrid( 1, 1 );

  if ( biasMul )
    cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageMulKernel<<<dimGrid,dimBlock>>>( output, input, biasMul, numberOfPixels );

  if ( biasAdd )
    cmtkEntropyMinimizationIntensityCorrectionFunctionalCUDAUpdateOutputImageAddKernel<<<dimGrid,dimBlock>>>( output, input, biasAdd, numberOfPixels );
}

