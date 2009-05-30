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

#ifndef __cmtkOverlapMeasures_h_included_
#define __cmtkOverlapMeasures_h_included_

#include <cmtkconfig.h>

#include <cmtkTypedArray.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{
/// Class for overlap measures between multiple segmentations.
class OverlapMeasures
{
public:
  /// Constructor: allocate local data structures and do precomputations (e.g., count labels, etc).
  OverlapMeasures( const std::vector<TypedArray::SmartPtr>& dataArrays );

  /// Enumeration for different region weightings.
  typedef enum
  {
    /// Equal weighting of all regions.
    Equal,
    /// Weighting of regions proportional to volume.
    Volume,
    /// Weighting of regions proportional to inverse volume.
    VolumeInverse
  } RegionWeightingMode;

  /** Compute groupwise overlap with advanced options.
   *\return Number of labels included in computation. If this is zero, the resulting overlap values are invalid.
   */
  size_t ComputeGroupwiseOverlap
  ( const int firstLabel, //!< Analysis starts with this label.
    const int numberOfLabels, //!< Analysis covers these labels
    double& overlapEqualWeighted, //!< Equal-weighted overlap score is returned herein.
    double& overlapVolumeWeighted, //!< Volume-weighted overlap score is returned herein. 
    double& overlapInverseWeighted //!< Inverse volume-weighted overlap score is returned herein.
    ) const;

  /** Compute simple groupwise overlap.
   *\return Number of labels included in computation. If this is zero, the resulting overlap values are invalid.
   */
  size_t ComputeGroupwiseOverlap
  ( double& overlapEqualWeighted, //!< Equal-weighted overlap score is returned herein.
    double& overlapVolumeWeighted, //!< Volume-weighted overlap score is returned herein. 
    double& overlapInverseWeighted //!< Inverse volume-weighted overlap score is returned herein.
    ) const
  {
    return this->ComputeGroupwiseOverlap( 0, this->m_MaxLabelValue+1, overlapEqualWeighted, overlapVolumeWeighted, overlapInverseWeighted );
  }

  /// Return maximum label value used in data.
  unsigned int GetMaxLabelValue() const
  {
    return this->m_MaxLabelValue;
  }

private:
  /// Number of images.
  size_t m_NumberOfImages;

  /// Number of pixels: the minimum number over all images.
  size_t m_NumberOfPixels;

  /// Maximum label value used in the data.
  unsigned int m_MaxLabelValue;

  /// Data arrays.
  std::vector<TypedArray::SmartPtr> m_DataArrays;

  /// Compute pairwise overlap minimum.
  double ComputePairwiseOverlapMinMax( double& overlap_min, double& overlap_max, const TypedArray::SmartPtr& data0, const TypedArray::SmartPtr& data1, const int label ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkOverlapMeasures_h_included_

