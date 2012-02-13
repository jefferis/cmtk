/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkLabelCombinationLocalBinaryShapeBasedAveraging_h_included_
#define __cmtkLabelCombinationLocalBinaryShapeBasedAveraging_h_included_

#include <cmtkconfig.h>

#include <Segmentation/cmtkLabelCombinationLocalWeighting.h>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Combination of binary segmentations by locally-weighted shape-based averaging.
 *\attention All labels maps are treated as binary maps, i.e., all labels not equal to zero are considered equal.
 */
class LabelCombinationLocalBinaryShapeBasedAveraging
  : public LabelCombinationLocalWeighting
{
public:
  /// This class.
  typedef LabelCombinationLocalBinaryShapeBasedAveraging Self;

  /// Parent class.
  typedef LabelCombinationLocalWeighting Superclass;

  /// Constructor: compute label combination.
  LabelCombinationLocalBinaryShapeBasedAveraging( const UniformVolume::SmartConstPtr targetImage ) : Superclass( targetImage ) {}

  /// Set flag to detect local outliers at each pixel in the co-registered distance maps.
  void SetDetectLocalOutliers( const bool detectOutliers = true )
  {
    this->m_DetectLocalOutliers = detectOutliers;
  }
  
  /// Add an atlas (pair of reformatted, target-matched intensity image and label map).
  void AddAtlas( const UniformVolume::SmartConstPtr image, const UniformVolume::SmartConstPtr atlas );

  /// Get resulting combined segmentation.
  virtual TypedArray::SmartPtr GetResult() const;  

private:
  /// Compute result for a region.
  void ComputeResultForRegion( const Self::TargetRegionType& region, TypedArray& result ) const;

  /// Signed distance maps for the atlas label maps.
  std::vector<UniformVolume::SmartConstPtr> m_AtlasDMaps;

  /// Flag for outlier detection.
  bool m_DetectLocalOutliers;

protected:
  /** Delete atlas with given index. 
   * Call inherited member, then delete distance map.
   */
  virtual void DeleteAtlas( const size_t i )
  {
    this->Superclass::DeleteAtlas( i );
    this->m_AtlasDMaps.erase( this->m_AtlasDMaps.begin() + i );
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationLocalBinaryShapeBasedAveraging_h_included_
