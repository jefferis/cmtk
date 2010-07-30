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

#include "cmtkSimpleLevelsetDevice_kernels.h"

#include "GPU/cmtkCUDA.h"
#include "GPU/cmtkDeviceMemory.h"
#include "GPU/cmtkSumReduction_kernel.h"

#include <cuda_runtime_api.h>

__global__
void
cmtkSimpleLevelsetDeviceUpdateInsideOutsideKernel( float* partialInsideSum, float* partialOutsideSum, int* partialInside, float* levelset, float* data, const int nPixels )
{
  int nInside = 0;
  float insideSum = 0;
  float outsideSum = 0;

  for ( int idx = threadIdx.x; idx < nPixels; idx += blockDim.x )
    {
      const float l = levelset[idx];
      const float d = data[idx];
      const int flag = (l>0) ? 1 : 0;

      nInside += flag;
      insideSum += flag*d;
      outsideSum += (1-flag)*d;
    }
  
  partialInside[threadIdx.x] = nInside;
  partialInsideSum[threadIdx.x] = insideSum;
  partialOutsideSum[threadIdx.x] = outsideSum;
}

void
cmtk::SimpleLevelsetDeviceUpdateInsideOutside( float* levelset, float* data, const int nPixels, float* insideSum, float* outsideSum, int* nInside )
{
  const int nThreads = 512;

  DeviceMemory<int> partialInside( nThreads );
  DeviceMemory<float> partialInsideSum( nThreads );
  DeviceMemory<float> partialOutsideSum( nThreads );
  
  cmtkSimpleLevelsetDeviceUpdateInsideOutsideKernel<<<1,nThreads>>>( partialInsideSum.Ptr(), partialOutsideSum.Ptr(), partialInside.Ptr(), levelset, data, nPixels );
  cmtkCheckLastErrorCUDA;
  
  *nInside = SumReduction( partialInside.Ptr(), nThreads );
  *insideSum = SumReduction( partialInsideSum.Ptr(), nThreads );
  *outsideSum = SumReduction( partialOutsideSum.Ptr(), nThreads );
}

__global__
void
cmtkSimpleLevelsetDeviceUpdateInsideOutsideKernel( float* levelset, float* data, const int nPixels, const float mInside, const float mOutside, const float ratioInOut, const float timeDelta, const float levelsetThreshold )
{
  for ( size_t n = threadIdx.x; n < nPixels; n += blockDim.x )
    {
      const float d = data[n];
      const float l = levelset[n];
      
      const float zInside = fabsf( mInside - d );
      const float zOutside = fabs( mOutside - d );
      
      const float delta = ( zInside>zOutside ) ? -timeDelta * ratioInOut : timeDelta / ratioInOut;
      levelset[n] = fminf( levelsetThreshold, fmaxf( -levelsetThreshold, l+delta ) );
    }
}

void
cmtk::SimpleLevelsetDeviceUpdateLevelset( float* levelset, float* data, const int nPixels, const float mInside, const float mOutside, const float ratioInOut, const float timeDelta, const float levelsetThreshold )
{
  cmtkSimpleLevelsetDeviceUpdateInsideOutsideKernel<<<1,512>>>( levelset, data, nPixels, mInside, mOutside, ratioInOut, timeDelta, levelsetThreshold );
  cmtkCheckLastErrorCUDA;
}
