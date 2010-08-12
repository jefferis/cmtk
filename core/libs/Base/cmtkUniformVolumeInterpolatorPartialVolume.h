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

#ifndef __cmtkUniformVolumeInterpolatorPartialVolume_h_included_
#define __cmtkUniformVolumeInterpolatorPartialVolume_h_included_

#include <cmtkconfig.h>

#include "Base/cmtkUniformVolumeInterpolatorBase.h"

#include "Base/cmtkVector3D.h"
#include "Base/cmtkUniformVolume.h"
#include "System/cmtkSmartPtr.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{
/** Partial-volume interpolation class.
 *
 */
class UniformVolumeInterpolatorPartialVolume :
  /// Inherit interface from base class.
  public UniformVolumeInterpolatorBase
{
public:
  /// This class type.
  typedef UniformVolumeInterpolatorPartialVolume Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  UniformVolumeInterpolatorPartialVolume( const UniformVolume& volume )
    : UniformVolumeInterpolatorBase( volume )
  {
  }
  
  /** Get data at location.
   *
   * This function performs interpolation of one value from m_Volume at location
   * v using the interpolation function given as the class template parameter.
   *
   * This function should return true if a value can be interpolated from
   * m_Volume at v, and it should return false if v is outside the range
   * where a value can be interpolated (i.e., outside the volume boundaries).
   */
  virtual bool GetDataAt( const Vector3D& v, Types::DataItem& value ) const;

  /// Get data at a pre-computed relative pixel index. This is faster if we already know the pixel index and fractional coordinate of a location.
  virtual Types::DataItem GetDataDirect( const int* imageGridPoint, const Types::Coordinate* insidePixel ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkUniformVolumeInterpolatorPartialVolume_h_included_

