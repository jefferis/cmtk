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

#include "Base/cmtkFilterMask.h"
#include "Base/cmtkUnits.h"

#include "System/cmtkProgress.h"

#include "GPU/cmtkDeviceMemory.h"
#include "GPU/cmtkDeviceUniformVolumeArray.h"
#include "GPU/cmtkDeviceImageConvolution_kernels.h"
#include "GPU/cmtkDeviceThresholdData_kernels.h"

void
cmtk::SimpleLevelsetDevice
::Evolve( const int numberOfIterations, const bool forceIterations )
{
//  const double normFactor = 1.0/(sqrt(2*M_PI) * sigma);
//  for ( size_t i = 0; i < radius; ++i )
//    {
//    this->m_HistogramKernel[idx][i] = cmtk::ScaleHistogramValueTrait<HistogramBinType>::Scale( normFactor * exp( -MathUtil::Square( 1.0 * i / sigma ) / 2 ) );
//    }

  const size_t numberOfPixels = this->m_Levelset->GetNumberOfPixels();

  DeviceUniformVolumeArray::SmartPtr deviceVolume = DeviceUniformVolumeArray::Create( *(this->m_Volume) );
  DeviceUniformVolumeArray::SmartPtr deviceLevelset = DeviceUniformVolumeArray::Create( *(this->m_Levelset) );
  
  DeviceMemory<float>::SmartPtr temporary = DeviceMemory<float>::Create( numberOfPixels );

  size_t nInsideOld = 0, nInside = 1;

  Progress::Begin( 0, numberOfIterations, 1, "Levelset Evolution" );
  for ( int it = 0; (it < numberOfIterations) && ((nInside!=nInsideOld) || forceIterations); ++it )
    {
    Progress::SetProgress( it );

//    cmtkDeviceImageConvolutionInPlace();

//    cmtkSimpleLevelsetDeviceStage1();
//    cmtkSimpleLevelsetDeviceStage2();

    cmtkDeviceThresholdData( temporary->Ptr(), numberOfPixels, -this->m_LevelsetThreshold, this->m_LevelsetThreshold );
    }

  temporary->CopyToHost( this->m_Levelset->GetData()->GetDataPtr(), numberOfPixels );

  Progress::Done();
}
