/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkImagePairAffineRegistrationFunctionalDevice.h"

#include <GPU/cmtkImagePairAffineRegistrationFunctionalDevice_kernels.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairAffineRegistrationFunctionalDevice::ImagePairAffineRegistrationFunctionalDevice
( UniformVolume::SmartConstPtr& fixedVolume, UniformVolume::SmartConstPtr& movingVolume ) 
  : ImagePairRegistrationFunctional( fixedVolume, movingVolume ),
    m_FixedVolumeOnDevice( DeviceUniformVolumeArray::Create( *(fixedVolume) ) ),
    m_MovingVolumeOnDevice( DeviceUniformVolumeArray::Create( *(movingVolume) ) )
{
}

ImagePairAffineRegistrationFunctionalDevice::ReturnType
ImagePairAffineRegistrationFunctionalDevice::Evaluate()
{
  const AffineXform::MatrixType xformMatrix;

  float matrix[4][4];
  for ( size_t j = 0; j < 4; ++j )
    {
    for ( size_t i = 0; i < 4; ++i )
      {
      matrix[j][i] = static_cast<float>( xformMatrix[j][i] );
      }
    }

  FixedVector<3,float> deltas = this->m_ReferenceGrid->Deltas();

  // multiply deltas for index-to-image space conversion
  for ( size_t j = 0; j < 3; ++j )
    {
    for ( size_t i = 0; i < 3; ++i )
      {
      matrix[j][i] *= deltas[j];
      }
    }

  // divide by size to get to normalized image coordinates after mirror
  for ( size_t j = 0; j < 4; ++j ) // here, need to run up to 3 because translation is also in output space
    {
    for ( size_t i = 0; i < 3; ++i )
      {
      matrix[j][i] /= static_cast<float>( this->m_FloatingGrid->m_Size[i] );
      }
    }
  
  return -ImagePairAffineRegistrationFunctionalDeviceEvaluateMSD( this->m_ReferenceGrid->m_Dims.begin(), this->m_FixedVolumeOnDevice->GetDeviceArrayPtr()->GetArrayOnDevice(), 
								  this->m_FloatingGrid->m_Dims.begin(), this->m_MovingVolumeOnDevice->GetDeviceArrayPtr()->GetArrayOnDevice(), matrix );
}

} // namespace cmtk
