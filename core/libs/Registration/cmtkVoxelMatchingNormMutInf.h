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

#ifndef __cmtkVoxelMatchingNormMutInf_h_included_
#define __cmtkVoxelMatchingNormMutInf_h_included_

#include <cmtkconfig.h>

#include "Registration/cmtkRegistrationJointHistogram.h"
#include "Base/cmtkInterpolator.h"
#include "System/cmtkSmartPtr.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Voxel metric "normalized mutual information".
 *\deprecated For future code, use cmtk::ImagePairSimilarityMetricNMI instead.
 */
template<Interpolators::InterpolationEnum I=Interpolators::LINEAR>
class VoxelMatchingNormMutInf : 
  /// Inherit basic functionality from 2D histogram.
  public RegistrationJointHistogram<I> 
{
public:
  /// This type.
  typedef VoxelMatchingNormMutInf<I> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;
  
  /// Parent class.
  typedef RegistrationJointHistogram<I> Superclass;

  /** Constructor.
   * For reference and floating volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param fltVolume The floating (transformed) volume.
   *@param numRefBins The desired number of bins to classify the 
   * reference data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   *@param numFltBins The desired number of bins to classify the 
   * floating image data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   */
  VoxelMatchingNormMutInf ( const UniformVolume* refVolume, const UniformVolume* fltVolume,
			    const unsigned int numRefBins = CMTK_HISTOGRAM_AUTOBINS, const unsigned int numFltBins = CMTK_HISTOGRAM_AUTOBINS )
    : RegistrationJointHistogram<I>( refVolume, fltVolume, numRefBins, numFltBins ) {};
  
  /** Constructor with explicit value range limits.
   * For reference and floating volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param fltVolume The floating (transformed) volume.
   *@param numRefBins The desired number of bins to classify the 
   * reference data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   *@param numFltBins The desired number of bins to classify the 
   * floating image data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   *@param minBoundRef Lower bound for reference image values; all actual pixel
   * values below this bound will be sorted into the lowest histogram bin.
   *@param minBoundRef Upper bound for reference image values; all actual pixel
   * values above this bound will be sorted into the highest histogram bin.
   *@param minBoundFlt Lower bound for floating image values; all actual pixel
   * values below this bound will be sorted into the lowest histogram bin.
   *@param minBoundFlt Upper bound for floating image values; all actual pixel
   * values above this bound will be sorted into the highest histogram bin.
   */
  VoxelMatchingNormMutInf ( const UniformVolume* refVolume, const UniformVolume* fltVolume, const Types::DataItemRange& rangeRef, const Types::DataItemRange& rangeFlt )
    : RegistrationJointHistogram<I>( refVolume, fltVolume, CMTK_HISTOGRAM_AUTOBINS, CMTK_HISTOGRAM_AUTOBINS, rangeRef, rangeFlt ) {};
  
  /// Return normalized mutual information.
  typename Self::ReturnType Get () const 
  {
    double HX, HY, HXY;
    
    this->GetMarginalEntropies(HX,HY);
    this->GetJointEntropy(HXY);
    
    return static_cast<typename Self::ReturnType>( (HX + HY) / HXY );
  }
};

/// Normalized mutual information with trilinear interpolation.
typedef VoxelMatchingNormMutInf<Interpolators::LINEAR> VoxelMatchingNormMutInf_Trilinear;

/// Normalized mutual information with nearest-neighbor interpolation.
typedef VoxelMatchingNormMutInf<Interpolators::NEAREST_NEIGHBOR> VoxelMatchingNormMutInf_NearestNeighbor;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVoxelMatchingNormMutInf_h_included_
