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

#ifndef __cmtkRegistrationJointHistogram_h_included_
#define __cmtkRegistrationJointHistogram_h_included_

#include <cmtkconfig.h>

#include <cmtkVoxelMatchingMetric.h>
#include <cmtkJointHistogram.h>
#include <cmtkInterpolator.h>

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
   */
  RegistrationJointHistogram
  ( const UniformVolume* refVolume, const UniformVolume* fltVolume,
    const unsigned int numBinsX = CMTK_HISTOGRAM_AUTOBINS, 
    const unsigned int numBinsY = CMTK_HISTOGRAM_AUTOBINS,
    const Types::DataItem minBoundX = -HUGE_VAL, const Types::DataItem maxBoundX = HUGE_VAL,
    const Types::DataItem minBoundY = -HUGE_VAL, const Types::DataItem maxBoundY = HUGE_VAL );
  
  /// Copy constructor.
  RegistrationJointHistogram ( RegistrationJointHistogram& other, const bool copyData = false );
  
  /// Copy constructor.
  RegistrationJointHistogram ( const RegistrationJointHistogram& other );

  unsigned int SetDataX( const UniformVolume* volume, const unsigned int numBins, const Types::DataItem minBound = -HUGE_VAL, const Types::DataItem maxBound = HUGE_VAL );

  unsigned int SetDataY( const UniformVolume* volume, const unsigned int numBins, const Types::DataItem minBound = -HUGE_VAL, const Types::DataItem maxBound = HUGE_VAL );
  
  void SetDataXY( const UniformVolume* volumeX, const unsigned int numBinsX, const UniformVolume* volumeY, const unsigned int numBinsY,
		  const Types::DataItem minBoundX = -HUGE_VAL, const Types::DataItem maxBoundX = HUGE_VAL, const Types::DataItem minBoundY = -HUGE_VAL, const Types::DataItem maxBoundY = HUGE_VAL );

  /** Continue incremental calculation by fractional voxel index.
   * For a given pair of reference and floating sample, the computation 
   * proceeds. This means the bin counters for both given samples are 
   * incremented. Important is that for increased efficiency, only the indices
   * of the samples in the original volume data are given. The bin indices of 
   * the respective data values are then retrieved from the pre-calculated raw
   * byte array (see InitDataset).
   *@param refIdx Index of the current reference data sample.
   *@param fltIdx Index of the current floating data sample.
   *@param location Fractional voxel coordinate of the probed floating data 
   * value.
   */
  void Proceed( const size_t refIdx, const size_t fltIdx, const Types::Coordinate* frac ) 
  {
#ifndef CMTK_PVI_HISTOGRAMS
    this->Increment( this->GetSampleX( refIdx ), this->GetSampleY( fltIdx, frac ) );
#else
    const byte refData = this->GetSampleX( refIdx );
    Types::Coordinate offsX = 1.0 - frac[0];
    Types::Coordinate offsY = 1.0 - frac[1];
    Types::Coordinate offsZ = 1.0 - frac[2];

    // set pointer to bottom-left-proximal node of enclosing voxel
    assert( (fltIdx + nextIJK) < NumSamplesY );
    byte *fltNode = DataY + fltIdx;

    // compute data by tri-linear interpolation
    this->Increment( refData, fltNode[0],       offsZ  * offsY  * offsX );
    this->Increment( refData, fltNode[1],       offsZ  * offsY  * frac[0] );
    this->Increment( refData, fltNode[nextJ],   offsZ  * frac[1]* offsX );
    this->Increment( refData, fltNode[nextIJ],  offsZ  * frac[1]* frac[0] );
    this->Increment( refData, fltNode[nextK],   frac[2]* offsY  * offsX );
    this->Increment( refData, fltNode[nextIK],  frac[2]* offsY  * frac[0] );
    this->Increment( refData, fltNode[nextJK],  frac[2]* frac[1]* offsX );
    this->Increment( refData, fltNode[nextIJK], frac[2]* frac[1]* frac[0] );
#endif
  }

  void AddMetric ( const Self& other )
  {
    this->AddJointHistogram( other );
  }

  void RemoveMetric ( const Self& other )
  {
    this->RemoveJointHistogram( other );
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkRegistrationJointHistogram_h_included_
