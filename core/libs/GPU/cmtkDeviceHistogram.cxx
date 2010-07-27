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

#include "cmtkDeviceHistogram.h"
#include "cmtkDeviceHistogram_kernels.h"

#include "System/cmtkException.h"
#include "System/cmtkMemory.h"

cmtk::DeviceHistogram
::DeviceHistogram( const size_t numberOfBins )
{
  this->m_NumberOfBins = numberOfBins;
  this->m_NumberOfBinsPadded = cmtk::Memory::GetNextPowerOfTwo( this->m_NumberOfBins );
  if ( this->m_NumberOfBinsPadded > 512 )
    {
    throw Exception( "Exceeded maximum number of histogram bins (512)" );
    }

  this->m_OnDeviceData = DeviceMemory<float>::Create( this->m_NumberOfBinsPadded );
  this->m_OnDeviceData->SetToZero();

  this->m_OnDeviceResult = DeviceMemory<float>::Create( this->m_NumberOfBinsPadded );
  this->m_OnDeviceResult->SetToZero();
}

void
cmtk::DeviceHistogram
::Reset()
{
  this->m_OnDeviceData->SetToZero();
}

void
cmtk::DeviceHistogram
::Populate( const DeviceMemory<float>& dataOnDevice, const float rangeFrom, const float rangeTo, const bool logScale )
{
  cmtkDeviceHistogramPopulate( this->m_OnDeviceData->Ptr(), dataOnDevice.Ptr(), rangeFrom, rangeTo, logScale, this->m_NumberOfBins, dataOnDevice.GetNumberOfItems() );
}

void
cmtk::DeviceHistogram
::Populate( const DeviceMemory<float>& dataOnDevice, const DeviceMemory<int>& maskOnDevice, const float rangeFrom, const float rangeTo, const bool logScale )
{
  cmtkDeviceHistogramPopulate( this->m_OnDeviceData->Ptr(), dataOnDevice.Ptr(), maskOnDevice.Ptr(), rangeFrom, rangeTo, logScale, this->m_NumberOfBins, dataOnDevice.GetNumberOfItems() );
}

float
cmtk::DeviceHistogram
::GetEntropy() const
{
  cmtkDeviceHistogramEntropy( this->m_OnDeviceResult->Ptr(), this->m_OnDeviceData->Ptr(), this->m_NumberOfBinsPadded );

  float result;
  this->m_OnDeviceResult->CopyToHost( &result, 1 );
  return result;
}
