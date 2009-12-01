/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#ifndef __cmtkUniformVolumeInterpolator_h_included_
#define __cmtkUniformVolumeInterpolator_h_included_

#include <cmtkconfig.h>

#include <cmtkUniformVolumeInterpolatorBase.h>

#include <cmtkSmartPtr.h>
#include <cmtkVector3D.h>
#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{
/** Class template for kernel-based volume interpolators.
 *
 * This class is templated over the interpolation function, e.g., linear, cubic, or sinc.
 *
 *\see LinearInterpolator
 *\see CubicInterpolator
 *\see SincInterpolator
 *\see NearestNeighborInterpolator
 */
template<class TInterpolationFunction>
class UniformVolumeInterpolator :
  /// Inherit interface from base class.
  public UniformVolumeInterpolatorBase
{
public:
  /// This class type.
  typedef UniformVolumeInterpolator<TInterpolationFunction> Self;

  /// Superclass type.
  typedef UniformVolumeInterpolatorBase Superclass;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  UniformVolumeInterpolator( const UniformVolume::SmartPtr& volume = UniformVolume::SmartPtr::Null )
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

  virtual Types::DataItem GetDataDirect( const size_t baseIndex, const int* imageGridPoint, const Types::Coordinate* insidePixel ) const;
};

//@}

} // namespace cmtk

#include <cmtkUniformVolumeInterpolator.txx>

#endif // #ifndef __cmtkUniformVolumeInterpolator_h_included_

