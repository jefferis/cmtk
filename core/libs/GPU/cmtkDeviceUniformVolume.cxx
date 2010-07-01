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

#include "cmtkDeviceUniformVolume.h"

cmtk::DeviceUniformVolume::
DeviceUniformVolume( const UniformVolume& volume, const size_t padDataToMultiple )
{
  // set volume parameters
  UniformVolumeOnDevice onDevice;
  for ( size_t i = 0; i < 3; ++i )
    {
    onDevice.m_Dims[i] = volume.m_Dims[i];
    onDevice.m_Delta[i] = volume.m_Delta[i];
    }

  // convert volume data to float and copy to device
  if ( volume.GetData() )
    {
    TypedArray::SmartPtr floatData = volume.GetData()->Convert( TYPE_FLOAT );

    this->m_OnDeviceData = DeviceMemory<float>::Create( volume.GetNumberOfPixels(), static_cast<float*>( floatData->GetDataPtr() ), padDataToMultiple );
    
    // set device pointer to data, then copy whole structure to device.
    onDevice.m_Data = this->m_OnDeviceData->Ptr();
    }

  this->m_OnDevice = DeviceMemory<UniformVolumeOnDevice>::Create( 1, &onDevice );
}
