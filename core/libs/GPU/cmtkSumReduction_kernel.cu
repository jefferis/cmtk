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

#include "cmtkSumReduction_kernel.h"

#include "GPU/cmtkCUDA.h"

#include <cuda_runtime_api.h>

template<class T>
__global__
void cmtkSumReductionKernel( T* data, const int n )
{
  const int tx = threadIdx.x;

  for ( int i = tx + blockDim.x; i < n; i += blockDim.x )
    {
      data[tx] += data[i];
    }

  __syncthreads();

  if ( tx == 0 )
    {
      for ( int i = 1; i < blockDim.x; ++i )
	data[0] += data[i];
    }
}

template<class T>
T
cmtk::SumReduction( T* data, const int n )
{
  cmtkSumReductionKernel<T><<<1,512>>>( data, n );
  
  T result;
  cmtkCheckCallCUDA( cudaMemcpy( &result, data, sizeof( T ), cudaMemcpyDeviceToHost ) );
  return result;
}

template int cmtk::SumReduction<int>( int* data, const int n );
template float cmtk::SumReduction<float>( float* data, const int n );
