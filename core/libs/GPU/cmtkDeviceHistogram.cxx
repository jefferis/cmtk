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

#include "cmtkException.h"

#include <limits.h>

cmtk::DeviceHistogram
::DeviceHistogram( const size_t numberOfBins )
{
  this->m_NumberOfBinsPadded = Self::GetNextPowerOfTwo( numberOfBins );
  if ( this->m_NumberOfBinsPadded > 512 )
    {
    throw Exception( "Exceeded maximum number of histogram bins (512)" );
    }

  this->m_OnDeviceData = DeviceMemory<float>::Create( this->m_NumberOfBinsPadded );
  this->m_OnDeviceData->SetToZero();

  this->m_OnDeviceResult = DeviceMemory<float>::Create( this->m_NumberOfBinsPadded );
  this->m_OnDeviceResult->SetToZero();
}

float
cmtk::DeviceHistogram
::GetEntropy() const
{
  cmtkDeviceHistogramEntropy( this->m_OnDeviceResult->Ptr(), this->m_OnDeviceData->Ptr(), this->m_NumberOfBinsPadded );

  float result;
  this->m_OnDeviceResult->CopyFromDevice( &result, 1 );
  return result;
}

size_t
cmtk::DeviceHistogram
::GetNextPowerOfTwo( size_t k )
{

// http://en.wikipedia.org/wiki/Power_of_two#Algorithm_to_find_the_next-highest_power_of_two 

  if (k == 0)
    return 1;
  
  k--;
  for (int i=1; i<sizeof(size_t)*CHAR_BIT; i<<=1)
    k = k | k >> i;

  return k+1;
}
