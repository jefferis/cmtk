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

#include "cmtkEntropyMinimizationIntensityCorrectionFunctionalBaseCUDA.h"

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalBaseCUDA
::SetInputImage( UniformVolume::SmartConstPtr& inputImage )
{
  this->Superclass::SetInputImage( inputImage );

  this->m_BiasFieldAddCUDA = DeviceMemoryCUDA<float>::Create( inputImage->GetNumberOfPixels() );
  this->m_BiasFieldMulCUDA = DeviceMemoryCUDA<float>::Create( inputImage->GetNumberOfPixels() );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalBaseCUDA
::SetBiasFieldAdd( const UniformVolume& biasFieldAdd )
{
  this->Superclass::SetBiasFieldAdd( biasFieldAdd );
  this->m_BiasFieldAddCUDA->CopyToDevice( this->m_BiasFieldAdd->GetDataPtrTemplate(), this->m_BiasFieldAdd->GetDataSize() );
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalBaseCUDA
::SetBiasFieldMul( const UniformVolume& biasFieldMul )
{
  this->Superclass::SetBiasFieldMul( biasFieldMul );
  this->m_BiasFieldMulCUDA->CopyToDevice( this->m_BiasFieldMul->GetDataPtrTemplate(), this->m_BiasFieldMul->GetDataSize() );
}
  
void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalBaseCUDA
::UpdateBiasFields( const bool foregroundOnly )
{
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalBaseCUDA
::UpdateBiasFieldAdd( const bool foregroundOnly )
{
}

void
cmtk::EntropyMinimizationIntensityCorrectionFunctionalBaseCUDA
::UpdateBiasFieldMul( const bool foregroundOnly )
{
}
