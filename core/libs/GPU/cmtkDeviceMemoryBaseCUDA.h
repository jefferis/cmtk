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

#ifndef __cmtkDeviceMemoryBaseCUDA_h_included_
#define __cmtkDeviceMemoryBaseCUDA_h_included_

#include <cmtkconfig.h>

#include <cmtkSmartConstPtr.h>
#include <cmtkSmartPtr.h>

#include <new>

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Resource managing class for raw memory allocated on a GPU device through CUDA.
class DeviceMemoryBaseCUDA
{
public:
  /// This class.
  typedef DeviceMemoryBaseCUDA Self;

  /// Smart pointer-to-const.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Exception for failed allocation.
  class bad_alloc : public std::bad_alloc {};
  
  /// Create new object and allocate memory.
  Self::SmartPtr Alloc( const size_t nItems, const size_t itemSize )
  {
    return Self::SmartPtr( new Self( nItems, itemSize ) );
  }
  
  /// Destructor: free memory through CUDA.
  ~DeviceMemoryBaseCUDA();

  /// Copy from host to device memory.
  void CopyToDevice( const void *const srcPtrHost, const size_t count );
  
  /// Copy from device to host memory.
  void CopyFromDevice( void *const dstPtrHost, const size_t count ) const;
  
  /// Copy between two device memory locations.
  void CopyToDevice( const Self& srcPtrDevice, const size_t count );
  
protected:
  /// Constructor: allocate memory through CUDA.
  DeviceMemoryBaseCUDA( const size_t n /**!< Number of items to allocate */, const size_t size /**!< Item size, e.g., sizeof( float ) */ );

  /** Raw pointer to allocated device memory.
   * Note that this is a device memory space pointer, which is not valid in
   * host memory and can, therefore, not be dereferenced in host code.
   */
  void* m_PointerDevice;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeviceMemoryBaseCUDA_h_included_
