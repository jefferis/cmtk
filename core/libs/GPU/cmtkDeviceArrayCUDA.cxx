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
#include <cuda_runtime.h>

#include "GPU/cmtkCUDA.h"

cmtk::DeviceArrayCUDA
::DeviceArrayCUDA( const FixedVector<3,int>& dims3 )
  : m_Dims( dims3 )
{
  const struct cudaChannelFormatDesc desc = cudaCreateChannelDesc<float>();
  cmtkCheckCallCUDA( cudaMalloc3DArray( &(this->m_DeviceArrayPtr), &desc, make_cudaExtent( this->m_Dims[0], this->m_Dims[1], this->m_Dims[2] ) ) );
}


cmtk::DeviceArrayCUDA
::~DeviceArrayCUDA()
{
  if ( this->m_DeviceArrayPtr )
    cmtkCheckCallCUDA( cudaFreeArray( this->m_DeviceArrayPtr ) );
}

void
cmtk::DeviceArrayCUDA
::CopyToDevice( const float* data )
{
  cudaMemcpy3DParms copyParams = {0};
  
  copyParams.srcPtr   = make_cudaPitchedPtr( (void*)data, this->m_Dims[0]*sizeof(float), this->m_Dims[0], this->m_Dims[1] );  
  copyParams.dstArray = this->m_DeviceArrayPtr;
  copyParams.extent   = make_cudaExtent( this->m_Dims[0], this->m_Dims[1], this->m_Dims[2] );
  copyParams.kind     = cudaMemcpyHostToDevice;

  cmtkCheckCallCUDA( cudaMemcpy3D( &copyParams ) );
}

void
cmtk::DeviceArrayCUDA
::CopyOnDevice( const float* data )
{
  cudaMemcpy3DParms copyParams = {0};
  
  copyParams.srcPtr   = make_cudaPitchedPtr( (void*)data, this->m_Dims[0]*sizeof(float), this->m_Dims[0], this->m_Dims[1] );  
  copyParams.dstArray = this->m_DeviceArrayPtr;
  copyParams.extent   = make_cudaExtent( this->m_Dims[0], this->m_Dims[1], this->m_Dims[2] );
  copyParams.kind     = cudaMemcpyDeviceToDevice;

  cmtkCheckCallCUDA( cudaMemcpy3D( &copyParams ) );
}
