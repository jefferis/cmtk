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

#ifndef __cmtkDeviceTransformedVolumeAxes_h_included_
#define __cmtkDeviceTransformedVolumeAxes_h_included_

#include <cmtkconfig.h>

#include "GPU/cmtkDeviceMemory.h"
#include "Base/cmtkTransformedVolumeAxes.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/// Device memory representation of a uniform volume with static helper functions.
class DeviceTransformedVolumeAxes
{
public:
  /// This class.
  typedef DeviceTransformedVolumeAxes Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Create device representation of volume object.
  static Self::SmartPtr Create( const TransformedVolumeAxes& axesOnHost )
  {
    return Self::SmartPtr( new Self( axesOnHost ) );
  }
  
private:
  /// Constructor.
  DeviceTransformedVolumeAxes( const TransformedVolumeAxes& axesOnHost );

  /// Managed device memory pointer to histogram data.
  FixedVector<3,DeviceMemory<float>::SmartPtr> m_OnDeviceData;
};

} // namespace cmtk

#endif // #ifndef __cmtkDeviceTransformedVolumeAxes_h_included_
