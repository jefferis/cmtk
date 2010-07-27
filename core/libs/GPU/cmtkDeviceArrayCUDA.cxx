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

#include "cmtkDeviceArrayCUDA.h"

#include <cuda_runtime_api.h>

#include <cstdio>

cmtk::DeviceArrayCUDA
::DeviceArrayCUDA( const FixedVector<3,int>& dims3 )
  : m_Dims( dims3 )
{
  const struct cudaChannelFormatDesc desc = cudaCreateChannelDesc( 32, 0, 0, 0, cudaChannelFormatKindFloat );
  
  cudaError_t cudaError = cudaGetLastError();
  if ( cudaError != cudaSuccess )
    {
    fprintf( stderr, "ERROR: cudaCreateChannelDesc() failed with error '%s'\n", cudaGetErrorString( cudaError ) );
    exit( 1 );      
    }
  
  struct cudaExtent extent;
  extent.width = this->m_Dims[0];
  extent.height = this->m_Dims[1];
  extent.depth = this->m_Dims[2];
  
  cudaError = cudaMalloc3DArray( &(this->m_DeviceArrayPtr), &desc, extent );
  if ( cudaError != cudaSuccess )
    {
    fprintf( stderr, "ERROR: cudaMalloc3DArray() failed with error '%s'\n", cudaGetErrorString( cudaError ) );
    this->m_DeviceArrayPtr = NULL;
    throw( Self::bad_alloc() );
    }
}


cmtk::DeviceArrayCUDA
::~DeviceArrayCUDA()
{
  if ( this->m_DeviceArrayPtr )
    cudaFreeArray( this->m_DeviceArrayPtr );
}

void
cmtk::DeviceArrayCUDA
::CopyToDevice( const float* data )
{
  cudaMemcpyToArray( this->m_DeviceArrayPtr, 0, 0, data, this->m_Dims.Sum(), cudaMemcpyHostToDevice);
}
