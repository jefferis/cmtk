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

#ifndef __cmtkRegistrationJointHistogram_h_included_
#define __cmtkRegistrationJointHistogram_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkJointHistogram.h>
#include <Base/cmtkInterpolator.h>

#include <Registration/cmtkVoxelMatchingMetric.h>

#ifndef CMTK_HISTOGRAM_AUTOBINS
/// Constant for number of bins to be determined automatically.
#define CMTK_HISTOGRAM_AUTOBINS 0
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{
/** 2-D histogram for entropy-based image similarity measures.
 */
template<Interpolators::InterpolationEnum I=Interpolators::LINEAR>
class RegistrationJointHistogram : 
  /// Inherit histogram with integral bins.
  public JointHistogram<int>,
  /// Inherit from generic voxel metric with internal byte data.
  public VoxelMatchingMetric<byte,TYPE_BYTE,I>
{
public:
  /// This class.
  typedef RegistrationJointHistogram<I> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Constructor.
   *@param refVolume The reference (fixed) volume.
   *@param fltVolume The floating (transformed) volume.
   *@param numBinsX The desired number of bins to classify the 
   * reference data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   *@param numBinsY The desired number of bins to classify the 
   * floating data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   *\param boundsX Value range for the X data distribution. Values outside this range will be assigned to the first and last histogram bins, respectively.
   *\param boundsY Value range for the Y data distribution. Values outside this range will be assigned to the first and last histogram bins, respectively.
   */
  RegistrationJointHistogram
  ( const UniformVolume* refVolume, const UniformVolume* fltVolume,
    const unsigned int numBinsX = CMTK_HISTOGRAM_AUTOBINS, 
    const unsigned int numBinsY = CMTK_HISTOGRAM_AUTOBINS,
    const Types::DataItemRange& boundsX = Types::DataItemRange( -HUGE_VAL, HUGE_VAL ),
    const Types::DataItemRange& boundsY = Types::DataItemRange( -HUGE_VAL, HUGE_VAL ) );
  
  /** Continue incremental calculation by fractional voxel index.
   * For a given pair of reference and floating sample, the computation 
   * proceeds. This means the bin counters for both given samples are 
   * incremented. Important is that for increased efficiency, only the indices
   * of the samples in the original volume data are given. The bin indices of 
   * the respective data values are then retrieved from the pre-calculated raw
   * byte array (see InitDataset).
   *@param refIdx Index of the current reference data sample.
   *@param fltIdx Index of the current floating data sample.
   *@param frac Fractional voxel coordinate of the probed floating data 
   * value.
   */
  inline void Proceed( const size_t refIdx, const size_t fltIdx, const Types::Coordinate* frac ) 
  {
    this->Increment( this->GetSampleX( refIdx ), this->GetSampleY( fltIdx, frac ) );
  }

  /// Add another metric object, e.g., from distributed computation.
  void AddMetric ( const Self& other )
  {
    this->AddJointHistogram( other );
  }

  /// Remove another metric, e.g., to undo results from an ROI in local computation.
  void RemoveMetric ( const Self& other )
  {
    this->RemoveJointHistogram( other );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRegistrationJointHistogram_h_included_
