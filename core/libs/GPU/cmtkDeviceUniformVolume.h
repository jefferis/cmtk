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

#ifndef __cmtkDeviceUniformVolume_h_included_
#define __cmtkDeviceUniformVolume_h_included_

#include <cmtkconfig.h>

#include "cmtkDeviceMemory.h"
#include "cmtkUniformVolumeOnDevice.h"

#include "Base/cmtkUniformVolume.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Device memory representation of a uniform volume with static helper functions.
class DeviceUniformVolume
{
public:
  /// This class.
  typedef DeviceUniformVolume Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Create device representation of volume object.
  static Self::SmartPtr Create( const UniformVolume& volume, const size_t padDataToMultiple = 1 /**!< Allocate device memory for data as multiple of this value.*/ )
  {
    return Self::SmartPtr( new Self( volume, padDataToMultiple ) );
  }
  
  /// Return device representation of volume.
  DeviceMemory<UniformVolumeOnDevice>& GetOnDevice()
  {
    return *(this->m_OnDevice);
  }

  /// Return device data pointer.
  DeviceMemory<float>& GetDataOnDevice()
  {
    return *(this->m_OnDeviceData);
  }

private:
  /// Constructor.
  DeviceUniformVolume( const UniformVolume& volume, const size_t padDataToMultiple = 1 );

  /// Managed device memory pointer to parameter block.
  DeviceMemory<UniformVolumeOnDevice>::SmartPtr m_OnDevice;

  /// Managed device memory pointer to volume data.
  DeviceMemory<float>::SmartPtr m_OnDeviceData;
};

} // namespace cmtk

#endif // #ifndef __cmtkDeviceUniformVolume_h_included_
