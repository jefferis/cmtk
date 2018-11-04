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

#include "cmtkDeviceMemoryCUDA.h"

#include <GPU/cmtkCUDA.h>

#include <cuda_runtime_api.h>

namespace cmtk {

DeviceMemoryCUDA ::DeviceMemoryCUDA(const size_t nBytes,
                                    const size_t padToMultiple) {
  this->m_NumberOfBytesAllocated =
      (((nBytes - 1) / padToMultiple) + 1) * padToMultiple;

  cmtkCheckCallCUDA(
      cudaMalloc(&(this->m_PointerDevice), this->m_NumberOfBytesAllocated));
}

DeviceMemoryCUDA ::~DeviceMemoryCUDA() {
  if (this->m_PointerDevice) cmtkCheckCallCUDA(cudaFree(this->m_PointerDevice));
}

void DeviceMemoryCUDA ::CopyToDevice(const void *const srcPtrHost,
                                     const size_t nBytes) {
  cmtkCheckCallCUDA(cudaMemcpy(this->m_PointerDevice, srcPtrHost, nBytes,
                               cudaMemcpyHostToDevice));
}

void DeviceMemoryCUDA ::CopyToHost(void *const dstPtrHost,
                                   const size_t nBytes) const {
  cmtkCheckCallCUDA(cudaMemcpy(dstPtrHost, this->m_PointerDevice, nBytes,
                               cudaMemcpyDeviceToHost));
}

void DeviceMemoryCUDA ::CopyOnDevice(const Self &srcPtrDevice,
                                     const size_t nBytes) {
  cmtkCheckCallCUDA(cudaMemcpy(this->m_PointerDevice,
                               srcPtrDevice.m_PointerDevice, nBytes,
                               cudaMemcpyDeviceToDevice));
}

void DeviceMemoryCUDA ::Memset(const int value, const size_t nBytes) {
  cmtkCheckCallCUDA(cudaMemset(this->m_PointerDevice, value, nBytes));
}

}  // namespace cmtk
