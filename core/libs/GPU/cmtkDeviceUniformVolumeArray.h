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

#ifndef __cmtkDeviceUniformVolumeTexture_h_included_
#define __cmtkDeviceUniformVolumeTexture_h_included_

#include <cmtkconfig.h>

#include "cmtkDeviceArray.h"

#include "Base/cmtkUniformVolume.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Representation of a uniform volume as 3D texture in device memory.
class DeviceUniformVolumeTexture
{
public:
  /// This class.
  typedef DeviceUniformVolumeTexture Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Create device representation of volume object.
  static Self::SmartPtr Create( const UniformVolume& volume )
  {
    return Self::SmartPtr( new Self( volume ) );
  }

  /// Get volume array on device.
  DeviceArray::DeviceArrayPointer GetDeviceArrayPtr()
  {
    return this->m_DeviceArrayPointer;
  }
  
private:
  /// Constructor.
  DeviceUniformVolumeTexture( const UniformVolume& volume );

  /// Device array pointer.
  DeviceArray::DeviceArrayPointer m_DeviceArrayPointer;
};

} // namespace cmtk

#endif // #ifndef __cmtkDeviceUniformVolumeTexture_h_included_
