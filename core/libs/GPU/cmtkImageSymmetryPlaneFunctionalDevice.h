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

#ifndef __cmtkImageSymmetryPlaneFunctionalDevice_h_included_
#define __cmtkImageSymmetryPlaneFunctionalDevice_h_included_

#include <cmtkconfig.h>

#include "Registration/cmtkImageSymmetryPlaneFunctionalBase.h"

#include "GPU/cmtkDeviceUniformVolumeArray.h"

namespace
cmtk
{

/** \addtogroup GPU */
//@{

/** Functional for finding a symmetry plane in 3-D volumes using GPU support.
 */
class ImageSymmetryPlaneFunctionalDevice :
  /// Inherit functional interface.
  public ImageSymmetryPlaneFunctionalBase
{
public:
  /// This class.
  typedef ImageSymmetryPlaneFunctionalDevice Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass
  typedef ImageSymmetryPlaneFunctionalBase Superclass;

  /// Constructor.
  ImageSymmetryPlaneFunctionalDevice( UniformVolume::SmartConstPtr& volume );

  /// Constructor with value range limits.
  ImageSymmetryPlaneFunctionalDevice( UniformVolume::SmartConstPtr& volume, const Types::DataItemRange& valueRange );

  /// Destructor.
  virtual ~ImageSymmetryPlaneFunctionalDevice() {}

  /// Compute functional value.
  virtual Self::ReturnType Evaluate();

private:
  /// Volume on compute device.
  DeviceUniformVolumeArray::SmartPtr m_VolumeOnDevice;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageSymmetryPlaneFunctionalDevice_h_included_
