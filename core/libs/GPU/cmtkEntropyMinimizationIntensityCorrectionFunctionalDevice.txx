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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice.h"
#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalDevice_kernels.h"

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>
::SetInputImage( UniformVolume::SmartConstPtr& inputImage )
{
  this->Superclass::SetInputImage( inputImage );
  this->m_InputImageDevice = DeviceUniformVolume::Create( *inputImage, 512 );
  this->m_NumberOfPixels = inputImage->GetNumberOfPixels();
}

template<unsigned int NOrderAdd,unsigned int NOrderMul>
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalDevice<NOrderAdd,NOrderMul>
::UpdateOutputImage( const bool foregroundOnly )
{
  if ( !this->m_OutputDataDevice )
    this->m_OutputDataDevice = DeviceMemory<float>::Create( this->m_NumberOfPixels, 512 );

  float* input = this->m_InputImageDevice->GetDataOnDevice().Ptr();
  float* output = this->m_OutputDataDevice->Ptr();

  if ( Self::PolynomialTypeMul::NumberOfMonomials )
    {
    cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( input, output, NOrderMul, 1 /*multiply*/ );
    input = output; // if additive bias also, apply to output of multiplicative stage
    }

  if ( Self::PolynomialTypeAdd::NumberOfMonomials )
    {
    cmtkEntropyMinimizationIntensityCorrectionFunctionalDeviceUpdateOutputImage( input, output, NOrderAdd, 0 /*multiply*/ );
    }
}
