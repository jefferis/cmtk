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

#ifndef __cmtkDeviceArray_h_included_
#define __cmtkDeviceArray_h_included_

#include <cmtkconfig.h>

#ifdef CMTK_USE_CUDA
#  include "cmtkDeviceArrayCUDA.h"
#else
#  include "cmtkDeviceArrayCL.h"
#endif

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Resource managing class template for type-specific memory allocated on a GPU device through .
template<class DeviceArrayGPU>
class DeviceArrayTemplate :
    /// Inherit privately from raw pointer base class.
    private DeviceArrayGPU
{
public:
  /// This class.
  typedef DeviceArrayTemplate<DeviceArrayGPU> Self;
  
  /// Smart pointer-to-const.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Base class.
  typedef DeviceArrayGPU Superclass;

  /// Device array pointer type forwarded.
  typedef typename Superclass::DeviceArrayPointer DeviceArrayPointer;

  /// Create new object and allocate memory.
  static typename Self::SmartPtr Create( const FixedVector<3,int>& dims3 /**!< Array dimensions */ )
  {
    return typename Self::SmartPtr( new Self( dims3 ) );
  }
  
  /// Create new object, allocate, and initialize memory.
  static typename Self::SmartPtr Create( const FixedVector<3,int>& dims3, /**!< Array dimensions */
					 const float* initFrom /**!< Initialize from this region in host memory.*/ )
  {
    Self* newObject = new Self( dims3 );
    newObject->CopyToDevice( initFrom );
    return typename Self::SmartPtr( newObject );
  }
  
private:
  /// Constructor: allocate memory on device through base class.
  DeviceArrayTemplate(  const FixedVector<3,int>& dims3 ) 
    : DeviceArrayGPU( dims3 )
  {}
};

#ifdef CMTK_USE_CUDA
typedef DeviceArrayTemplate<DeviceArrayCUDA> DeviceArray;
#else
typedef DeviceArrayTemplate<DeviceArrayCL> DeviceArray;
#endif

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeviceArray_h_included_
