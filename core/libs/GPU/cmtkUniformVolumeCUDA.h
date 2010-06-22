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

#ifndef __cmtkUniformVolumeCUDA_h_included_
#define __cmtkUniformVolumeCUDA_h_included_

#include <cmtkconfig.h>

#include "cmtkDeviceMemoryCUDA.h"
#include "cmtkUniformVolumeOnDeviceCUDA.h"

#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Device memory representation of a uniform volume with static helper functions.
class UniformVolumeCUDA
{
public:
  /// This class.
  typedef UniformVolumeCUDA Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Create device representation of volume object.
  static Self::SmartPtr Create( const UniformVolume& volume )
  {
    return Self::SmartPtr( new Self( volume ) );
  }
  
  /// Return device representation of volume.
  DeviceMemoryCUDA<UniformVolumeOnDeviceCUDA>& GetOnDevice()
  {
    return *(this->m_OnDevice);
  }

  /// Return device data pointer.
  DeviceMemoryCUDA<float>& GetDataOnDevice()
  {
    return *(this->m_OnDeviceData);
  }

private:
  /// Constructor.
  UniformVolumeCUDA( const UniformVolume& volume );

  /// Managed device memory pointer to parameter block.
  DeviceMemoryCUDA<UniformVolumeOnDeviceCUDA>::SmartPtr m_OnDevice;

  /// Managed device memory pointer to volume data.
  DeviceMemoryCUDA<float>::SmartPtr m_OnDeviceData;
};

} // namespace cmtk

#endif // #ifndef __cmtkUniformVolumeCUDA_h_included_
