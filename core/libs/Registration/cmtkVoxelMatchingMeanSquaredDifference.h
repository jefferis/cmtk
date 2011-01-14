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

#ifndef __cmtkVoxelMatchingMeanSquaredDifference_h_included_
#define __cmtkVoxelMatchingMeanSquaredDifference_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkVoxelMatchingMetric.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkMathUtil.h>

#include <System/cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/** Mean squared difference metric.
 *\deprecated For future code, use cmtk::ImagePairSimilarityMetricMSD instead.
 */
class VoxelMatchingMeanSquaredDifference :
  /// Inherit generic voxel metric with internal short data.
  public VoxelMatchingMetricShort
{
public:
  /// This type.
  typedef VoxelMatchingMeanSquaredDifference Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Constructor.
   * For reference and model volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param fltVolume The floating (moving) volume.
   */
  VoxelMatchingMeanSquaredDifference( const UniformVolume* refVolume, const UniformVolume* fltVolume );

  /** Add a pair of values to the metric.
   */
  template<class T> void Increment( const T a, const T b )
  {
    if ( (a == this->DataX.padding()) || (b == this->DataY.padding()) ) return;
    ++Samples;
    Sum -= MathUtil::Square( a - b );
  }

  /** Remove a pair of values from the metric.
   */
  template<class T> void Decrement( const T a, const T b )
  {
    if ( (a == this->DataX.padding()) || (b == this->DataY.padding()) ) return;
    --Samples;
    Sum += MathUtil::Square( a - b );
  }

  /// Reset internal variables for next computation.
  void Reset () 
  {
    Sum = 0;
    Samples = 0;
  }

  /// Get similarity measure value.
  Self::ReturnType Get() const 
  {
    return static_cast<Self::ReturnType>( Sum / Samples );
  }

  void AddMetric ( const Self& other )
  {
    Sum += other.Sum;
    Samples += other.Samples;
  }

  void RemoveMetric ( const Self& other )
  {
    Sum -= other.Sum;
    assert( Sum <= 0 );
    Samples -= other.Samples;
    assert( Samples >= 0 );
  }

private:
  /// Sum of all value pair differences.
  double Sum;

  /// Counter for number of sample pairs.
  int Samples;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVoxelMatchingMeanSquaredDifference_h_included_
