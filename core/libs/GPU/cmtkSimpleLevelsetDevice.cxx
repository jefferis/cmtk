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

#include "cmtkSimpleLevelsetDevice.h"

#include "Base/cmtkGaussianKernel.h"
#include "Base/cmtkUnits.h"

#include "System/cmtkProgress.h"

#include "GPU/cmtkDeviceMemory.h"
#include "GPU/cmtkDeviceUniformVolume.h"
#include "GPU/cmtkDeviceUniformVolumeArray.h"
#include "GPU/cmtkDeviceImageConvolution_kernels.h"
#include "GPU/cmtkDeviceThresholdData_kernels.h"
#include "GPU/cmtkSimpleLevelsetDevice_kernels.h"

#include <vector>

void
cmtk::SimpleLevelsetDevice
::Evolve( const int numberOfIterations, const bool forceIterations )
{
  FixedVector< 3, std::vector<float> > kernels;

  for ( int dim = 0; dim < 3; ++dim )
    {
    kernels[dim] = GaussianKernel<float>::GetSymmetricKernel( this->m_FilterSigma / this->m_Volume->Deltas()[dim], 0.01 /*maxError*/ );
    }

  const size_t numberOfPixels = this->m_Levelset->GetNumberOfPixels();

  DeviceUniformVolume::SmartPtr deviceVolume = DeviceUniformVolume::Create( *(this->m_Volume) );
  DeviceUniformVolumeArray::SmartPtr deviceLevelset = DeviceUniformVolumeArray::Create( *(this->m_Levelset) );
  
  DeviceMemory<float>::SmartPtr temporary = DeviceMemory<float>::Create( numberOfPixels );

  int nInsideOld = 0, nInside = 1;

  Progress::Begin( 0, numberOfIterations, 1, "Levelset Evolution" );
  for ( int it = 0; (it < numberOfIterations) && ((nInside!=nInsideOld) || forceIterations); ++it )
    {
    Progress::SetProgress( it );

    DeviceImageConvolution( temporary->Ptr(), this->m_Volume->GetDims().begin(), deviceLevelset->GetDeviceArrayPtr()->GetArrayOnDevice(), 
			    kernels[0].size(), &kernels[0][0], kernels[1].size(), &kernels[1][0], kernels[2].size(), &kernels[2][0] );
    
    float insideSum, outsideSum;
    SimpleLevelsetDeviceUpdateInsideOutside( temporary->Ptr(), deviceVolume->GetDataOnDevice().Ptr(), numberOfPixels, &insideSum, &outsideSum, &nInside );

    const int nOutside = numberOfPixels - nInside;
    SimpleLevelsetDeviceUpdateLevelset( temporary->Ptr(), deviceVolume->GetDataOnDevice().Ptr(), numberOfPixels, insideSum / nInside, outsideSum / nOutside, 1.0 * nInside / nOutside, this->m_TimeDelta, this->m_LevelsetThreshold );
    }

  temporary->CopyToHost( this->m_Levelset->GetData()->GetDataPtr(), numberOfPixels );

  Progress::Done();
}
