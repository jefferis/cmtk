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

#ifndef __cmtkDeviceMemoryCUDA_h_included_
#define __cmtkDeviceMemoryCUDA_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>
#include <System/cmtkSmartConstPtr.h>
#include <System/cmtkSmartPtr.h>

#include <new>

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Resource managing class for raw memory allocated on a GPU device through CUDA.
class DeviceMemoryCUDA
    /// Make sure this is never copied.
  : private CannotBeCopied
{
public:
  /// This class.
  typedef DeviceMemoryCUDA Self;

  /// Smart pointer-to-const.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Exception for failed allocation.
  class bad_alloc : public std::bad_alloc {};
  
  /// Constructor: allocate memory through CUDA.
  DeviceMemoryCUDA( const size_t nBytes /*!< Number of bytes to allocate */, const size_t padToMultiple = 1 /*!< Pad to allocate nearest multiple of this many bytes. */ );

  /// Destructor: free memory through CUDA.
  virtual ~DeviceMemoryCUDA();

  /// Get number of bytes allocated on device.
  size_t GetNumberOfBytesAllocated()
  {
    return this->m_NumberOfBytesAllocated;
  }

protected:
  /// Create new object and allocate memory.
  Self::SmartPtr Alloc( const size_t nBytes, const size_t padToMultiple = 1 )
  {
    return Self::SmartPtr( new Self( nBytes, padToMultiple ) );
  }
  
  /// Copy from host to device memory.
  void CopyToDevice( const void *const srcPtrHost, const size_t nBytes );
  
  /// Copy from device to host memory.
  void CopyToHost( void *const dstPtrHost, const size_t nBytes ) const;
  
  /// Copy between two device memory locations.
  void CopyOnDevice( const Self& srcPtrDevice, const size_t nBytes );
  
  /// Copy between two device memory locations.
  void Memset( const int value, const size_t nBytes );
  
  /** Raw pointer to allocated device memory.
   * Note that this is a device memory space pointer, which is not valid in
   * host memory and can, therefore, not be dereferenced in host code.
   */
  void* m_PointerDevice;

  /// Total number of bytes allocated on device.
  size_t m_NumberOfBytesAllocated;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeviceMemoryCUDA_h_included_
