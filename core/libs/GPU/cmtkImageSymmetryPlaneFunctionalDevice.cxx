/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkImageSymmetryPlaneFunctionalDevice.h"

#include "Base/cmtkTransformedVolumeAxes.h"

#include "GPU/cmtkImageSymmetryPlaneFunctionalDevice_kernels.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImageSymmetryPlaneFunctionalDevice::ImageSymmetryPlaneFunctionalDevice
( UniformVolume::SmartConstPtr& volume ) 
  : ImageSymmetryPlaneFunctionalBase( volume ),
    m_VolumeOnDevice( DeviceUniformVolumeTexture::Create( *(this->m_Volume) ) )
{
}

ImageSymmetryPlaneFunctionalDevice::ImageSymmetryPlaneFunctionalDevice
( UniformVolume::SmartConstPtr& volume, 
  const Types::DataItemRange& valueRange )
  : ImageSymmetryPlaneFunctionalBase( volume, valueRange ),
    m_VolumeOnDevice( DeviceUniformVolumeTexture::Create( *(this->m_Volume) ) )
{
}

ImageSymmetryPlaneFunctionalDevice::ReturnType
ImageSymmetryPlaneFunctionalDevice::Evaluate()
{
  return cmtkImageSymmetryPlaneFunctionalDeviceEvaluate( this->m_VolumeOnDevice->GetDeviceArrayPtr() );
}

} // namespace cmtk
