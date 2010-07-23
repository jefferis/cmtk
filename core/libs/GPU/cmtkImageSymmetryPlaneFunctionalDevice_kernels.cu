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

#include <cuda_runtime_api.h>

#include <cstdio>

/// Texture reference to volume data.
texture<float, 3, cudaReadModeElementType> texRef;

__constant__ float deviceAxesTNL[16384];

float
cmtkImageSymmetryPlaneFunctionalDeviceEvaluate( const int* dims3, void* array, const float* axesTNL )
{
  if ( (cudaMemcpyToSymbol( deviceAxesTNL, axesTNL, 3*(dims3[0]+dims3[1]+dims3[2])*sizeof( *axesTNL ), 0, cudaMemcpyHostToDevice ) != cudaSuccess) )
    {
      fprintf( stderr, "ERROR: cudaMemcpy() to constant memory failed with error %s\n",cudaGetErrorString( cudaGetLastError() ) );
      exit( 1 );      
    }
  
  // Set texture parameters
  texRef.addressMode[0] = cudaAddressModeWrap;
  texRef.addressMode[1] = cudaAddressModeWrap;
  texRef.addressMode[2] = cudaAddressModeWrap;
  texRef.filterMode = cudaFilterModeLinear; 
  texRef.normalized = true; 

  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  
  // Bind the array to the texture reference 
  cudaBindTextureToArray( texRef, (struct cudaArray*) array, channelDesc );

  cudaUnbindTexture( texRef );

  return 0;
}
