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

#ifndef __cmtkDeviceMemory_h_included_
#define __cmtkDeviceMemory_h_included_

#include <cmtkconfig.h>

#ifdef CMTK_USE_CUDA
#  include "cmtkDeviceMemoryCUDA.h"
#else
#  include "cmtkDeviceMemoryCL.h"
#endif

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Resource managing class template for type-specific memory allocated on a GPU device through .
#ifdef CMTK_USE_CUDA
template<typename T,class DeviceMemoryGPU = DeviceMemoryCUDA>
#else
template<typename T,class DeviceMemoryGPU = DeviceMemoryCL>
#endif
class DeviceMemory :
    /// Inherit privately from raw pointer base class.
    private DeviceMemoryGPU
{
public:
  /// This class.
  typedef DeviceMemory<T,DeviceMemoryGPU> Self;
  
  /// Smart pointer-to-const.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Base class.
  typedef DeviceMemoryGPU Superclass;

  /// Constructor: allocate memory on device through base class.
  DeviceMemory( const size_t n, /**!< Number of items.*/ const size_t padToMultiple = 1 ) 
    : DeviceMemoryGPU( n * sizeof( T ), padToMultiple * sizeof( T ) ),
      m_NumberOfItems( n )
  {}

  /// Create new object and allocate memory.
  static typename Self::SmartPtr Create( const size_t nItems, /**!< Allocate (at least) this many items of type T.*/ 
					 const size_t padToMultiple = 1 /**!< Pad number of allocated elements to next multiple of this number.*/  )
  {
    return typename Self::SmartPtr( new Self( nItems, padToMultiple ) );
  }
  
  /// Create new object, allocate, and initialize memory.
  static typename Self::SmartPtr Create( const size_t nItems, /**!< Allocate (at least) this many items of type T.*/ 
					 const T* initFrom, /**!< Initialize from this region in host memory.*/
					 const size_t padToMultiple = 1 /**!< Pad number of allocated elements to next multiple of this number.*/  )
  {
    Self* newObject = new Self( nItems, padToMultiple );
    newObject->CopyToDevice( initFrom, nItems );
    return typename Self::SmartPtr( newObject );
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
  void CopyToHost( void *const dstPtrHost, const size_t count ) const
  {
    this->Superclass::CopyToHost( dstPtrHost, count * sizeof( T ) );
  }
  
  /// Copy between two device memory locations.
  void CopyOnDevice( const Self& srcPtrDevice, const size_t count )
  {
    this->Superclass::CopyOnDevice( srcPtrDevice, count * sizeof( T ) );
  }

  /// Clear device memory (set to all zeroes).
  void SetToZero()
  {
    this->Superclass::Memset( 0, this->m_NumberOfItems * sizeof( T ) );
  }

  /// Get number of items.
  size_t GetNumberOfItems() const
  {
    return this->m_NumberOfItems;
  }

private:
  /// Number of items allocated.
  size_t m_NumberOfItems;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeviceMemory_h_included_
