/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkVoxelMatchingMetric_h_included_
#define __cmtkVoxelMatchingMetric_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkVoxelMatchingMetric_Type.h>

#include <System/cmtkSmartPtr.h>

#include <Base/cmtkInterpolator.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypes.h>
#include <Base/cmtkFunctional.h>

#include <cassert>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for voxel metrics with pre-converted image data.
 */
template<class T, ScalarDataType DT, Interpolators::InterpolationEnum I=cmtk::Interpolators::LINEAR>
class VoxelMatchingMetric :
  /// Inherit from type-template class.
  public VoxelMatchingMetric_Type<T,DT>
{
public:
  /// This type.
  typedef VoxelMatchingMetric<T,DT,I> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Return type: same as cmtk::Functional.
  typedef Functional::ReturnType ReturnType;

  /** Constructor.
   * For reference and model volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param fltVolume The model (transformed) volume.
   *@param initData If this flag is set (default), then the internal 
   * representation of the pixel data for both volumes is also created.
   * Derived classes may want to prevent this if they define their own
   * specific initialization (e.g., igsVoxelMatchingJointHistogram).
   */
  VoxelMatchingMetric( const UniformVolume* refVolume, const UniformVolume* fltVolume, const bool initData = true );

  /** Default constructor.
   */
  VoxelMatchingMetric() {};

  /// Get a value from the X distribution (reference image).
  T GetSampleX ( const size_t index ) const
  { 
    return this->DataX.Data[index]; 
  }
  
  /// Get a value from the Y distribution (floating image).
  T GetSampleY ( const size_t index ) const 
  { 
    return this->DataY.Data[index]; 
  }
  
  /// Interpolate a value from the Y distribution (floating image).
  T GetSampleY ( const size_t baseIndex, const Types::Coordinate* frac ) const;
};

/// Convenience typedef.
typedef VoxelMatchingMetric<short,TYPE_SHORT> VoxelMatchingMetricShort;

/// Convenience typedef.
typedef VoxelMatchingMetric<byte,TYPE_BYTE> VoxelMatchingMetricByte;

/// Convenience typedef.
typedef VoxelMatchingMetric<short,TYPE_SHORT,Interpolators::NEAREST_NEIGHBOR> VoxelMatchingMetricShort_NN;

/// Convenience typedef.
typedef VoxelMatchingMetric<byte,TYPE_BYTE,Interpolators::NEAREST_NEIGHBOR> VoxelMatchingMetricByte_NN;

//@}

} // namespace cmtk

#include "cmtkVoxelMatchingMetric.txx"

#endif // #ifndef __cmtkVoxelMatchingMetric_h_included_
