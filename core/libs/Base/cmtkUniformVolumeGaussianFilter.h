/*
//
//  Copyright 2004-2012 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkUniformVolumeGaussianFilter_h_included_
#define __cmtkUniformVolumeGaussianFilter_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkDataGridFilter.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkUnits.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Filter operations for 3D uniform image data.
 */
class UniformVolumeGaussianFilter :
  /// Prevent copying by inheritance.
  private DataGridFilter
{
public:
  /// This class.
  typedef UniformVolumeGaussianFilter Self;

  /// Constructor: link to UniformVolume object.
  explicit UniformVolumeGaussianFilter( UniformVolume::SmartConstPtr& volume ) : DataGridFilter( volume ), m_UniformVolume( volume ) {}

  /// Gaussian filter (using faster, separable filtering).
  TypedArray::SmartPtr GetFiltered3D( const Units::GaussianSigma& sigma, /*!< Kernel parameter "sigma" (standard deviation) */
				      const Types::Coordinate maxError = 0.01 /*!< Maximum approximation error: the kernel is truncated when it falls below this threshold */ ) const;
  
  /// Gaussian 1D filter.
  TypedArray::SmartPtr GetFiltered1D( const int direction /*!< Coordinate direction: 0 = x, 1 = y, 2 = z */,
				      const Units::GaussianSigma& sigma, /*!< Kernel parameter "sigma" (standard deviation) */
				      const Types::Coordinate maxError = 0.01 /*!< Maximum approximation error: the kernel is truncated when it falls below this threshold */ ) const;
  
private:
  /// The UniformVolume object we're working on.
  UniformVolume::SmartConstPtr m_UniformVolume;
};

} // namespace cmtk

#endif // #ifndef __cmtkUniformVolumeGaussianFilter_h_included_
