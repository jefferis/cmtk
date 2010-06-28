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

#ifndef __cmtkDeviceHistogram_h_included_
#define __cmtkDeviceHistogram_h_included_

#include <cmtkconfig.h>

#include "cmtkDeviceMemory.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Device memory representation of a uniform volume with static helper functions.
class DeviceHistogram
{
public:
  /// This class.
  typedef DeviceHistogram Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Create device representation of volume object.
  static Self::SmartPtr Create( const size_t numberOfBins /**!< Allocate device memory for data as multiple of this value.*/ )
  {
    return Self::SmartPtr( new Self( numberOfBins ) );
  }
  
  /// Return device data pointer.
  DeviceMemory<float>& GetDataOnDevice()
  {
    return *(this->m_OnDeviceData);
  }

  /// Get entropy.
  float GetEntropy() const;

private:
  /// Constructor.
  DeviceHistogram( const size_t numberOfBins );

  /// Number of bins after padding to power of 2.
  size_t m_NumberOfBinsPadded;

  /// Managed device memory pointer to histogram data.
  DeviceMemory<float>::SmartPtr m_OnDeviceData;

  /// Managed device memory pointer to result of histogram operations.
  mutable DeviceMemory<float>::SmartPtr m_OnDeviceResult;

  /** Utility function: get next power of two.
   *\url http://en.wikipedia.org/wiki/Power_of_two#Algorithm_to_find_the_next-highest_power_of_two 
   */
  static size_t GetNextPowerOfTwo( size_t k );
};

} // namespace cmtk

#endif // #ifndef __cmtkDeviceHistogram_h_included_
