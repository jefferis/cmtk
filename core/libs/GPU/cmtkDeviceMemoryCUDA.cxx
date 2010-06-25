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

#include "cmtkDeviceMemoryBaseCUDA.h"

#include <cuda_runtime_api.h>

namespace 
cmtk
{

DeviceMemoryBaseCUDA
::DeviceMemoryBaseCUDA( const size_t nBytes, const size_t padToMultiple )
{
  const size_t totalBytes = (((nBytes-1) / padToMultiple)+1) * padToMultiple;
  if ( cudaMalloc( &(this->m_PointerDevice), totalBytes ) != cudaSuccess )
    {
    this->m_PointerDevice = NULL;
    throw( Self::bad_alloc() );
    }
}

DeviceMemoryBaseCUDA
::~DeviceMemoryBaseCUDA()
{
  if ( this->m_PointerDevice )
    cudaFree( this->m_PointerDevice );
}

void
DeviceMemoryBaseCUDA
::CopyToDevice( const void *const srcPtrHost, const size_t nBytes )
{
  cudaMemcpy( this->m_PointerDevice, srcPtrHost, nBytes, cudaMemcpyHostToDevice );
}
  
void
DeviceMemoryBaseCUDA
::CopyFromDevice( void *const dstPtrHost, const size_t nBytes ) const
{
  cudaMemcpy( dstPtrHost, this->m_PointerDevice, nBytes, cudaMemcpyDeviceToHost );
} 

void
DeviceMemoryBaseCUDA
::CopyOnDevice( const Self& srcPtrDevice, const size_t nBytes )
{
  cudaMemcpy( this->m_PointerDevice, srcPtrDevice.m_PointerDevice, nBytes, cudaMemcpyDeviceToDevice );
}

} // namespace cmtk
