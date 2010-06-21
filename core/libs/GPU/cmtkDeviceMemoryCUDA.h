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

#include "cmtkDeviceMemoryBaseCUDA.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Resource managing class template for type-specific memory allocated on a GPU device through CUDA.
template<typename T>
class DeviceMemoryCUDA :
    /// Inherit privately from raw pointer base class.
    private DeviceMemoryBaseCUDA
{
public:
  /// This class.
  typedef DeviceMemoryCUDA<T> Self;
  
  /// Smart pointer-to-const.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Base class.
  typedef DeviceMemoryBaseCUDA Superclass;

  /// Create new object and allocate memory.
  static Self::SmartPtr Create( const size_t nItems )
  {
    return Self::SmartPtr( new Self( nItems ) );
  }
  
  /// Get const pointer.
  const T* Ptr() const
  {
    return static_cast<const T*>( this->m_PointerDevice );
  }

  /// Get non-const pointer.
  T* Ptr()
  {
    return static_cast<T*>( this->m_PointerDevice );
  }

  /// Copy from host to device memory.
  void CopyToDevice( const T *const srcPtrHost, const size_t count )
  {
    this->Superclass::CopyToDevice( srcPtrHost, count * sizeof( T ) );
  }
  
  /// Copy from device to host memory.
  void CopyFromDevice( void *const dstPtrHost, const size_t count ) const
  {
    this->Superclass::CopyFromDevice( dstPtrHost, count * sizeof( T ) );
  }
  
  /// Copy between two device memory locations.
  void CopyToDevice( const Self& srcPtrDevice, const size_t count )
  {
    this->Superclass::CopyToDevice( srcPtrDevice, count * sizeof( T ) );
  }
  
private:
  /// Constructor: allocate memory through CUDA.
  DeviceMemoryCUDA( const size_t n /**!< Number of items. */ ) : DeviceMemoryBaseCUDA( n, sizeof( T ) ) {}
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeviceMemoryCUDA_h_included_
