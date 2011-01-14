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

#ifndef __cmtkVoxelMatchingMutInf_h_included_
#define __cmtkVoxelMatchingMutInf_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkRegistrationJointHistogram.h>
#include <Base/cmtkInterpolator.h>
#include <System/cmtkSmartPtr.h>

#include <stdio.h>
#include <cassert>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Voxel metric "mutual information".
 *\deprecated For future code, use cmtk::ImagePairSimilarityMetricNCC instead.
 */
template<Interpolators::InterpolationEnum I=Interpolators::LINEAR>
class VoxelMatchingMutInf : 
  /// Inherit basic functionality from 2D histogram.
  public RegistrationJointHistogram<I> 
{
public:
  /// This type.
  typedef VoxelMatchingMutInf<I> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Constructor.
   * For reference and model volume, InitDataset is called.
   *@param refVolume The reference (fixed) volume.
   *@param fltVolume The model (transformed) volume.
   *@param numRefBins The desired number of bins to classify the 
   * reference data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   *@param numFltBins The desired number of bins to classify the 
   * model data. If this parameter is zero (default), a suitable value
   * is automatically determined.
   */
  VoxelMatchingMutInf ( const UniformVolume* refVolume, const UniformVolume* fltVolume, const unsigned int numRefBins = CMTK_HISTOGRAM_AUTOBINS, const unsigned int numFltBins = CMTK_HISTOGRAM_AUTOBINS )
    : RegistrationJointHistogram<I>( refVolume, fltVolume, numRefBins, numFltBins ) {};
  
  /// Return mutual information.
  typename Self::ReturnType Get () const 
  {
    double HX, HY;    
    this->GetMarginalEntropies(HX,HY);
    const double HXY = this->GetJointEntropy();
    
    return static_cast<typename Self::ReturnType>( HX + HY - HXY );
  }
};

/// Mutual information with trilinear interpolation.
typedef VoxelMatchingMutInf<Interpolators::LINEAR> VoxelMatchingMutInf_Trilinear;

/// Mutual information with nearest-neighbor interpolation.
typedef VoxelMatchingMutInf<Interpolators::NEAREST_NEIGHBOR> VoxelMatchingMutInf_NearestNeighbor;

//@}

} // namespace cmtk

#endif // #ifndef _VoxelMatchingMUTINF_H_
