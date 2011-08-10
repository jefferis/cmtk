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

#ifndef __cmtkLabelCombinationLocalShapeBasedAveraging_h_included_
#define __cmtkLabelCombinationLocalShapeBasedAveraging_h_included_

#include <cmtkconfig.h>

#include "cmtkLabelCombinationLocalVoting.h"

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Segmentation combination by locally-weighted shape-based averaging.
 *\attention Currently all labels maps are treated as binary maps, i.e., all labels not equal to zero are considered equal.
 */
class LabelCombinationLocalShapeBasedAveraging
  : public LabelCombinationLocalVoting
{
public:
  /// This class.
  typedef LabelCombinationLocalShapeBasedAveraging Self;

  /// Parent class.
  typedef LabelCombinationLocalVoting Superclass;

  /// Constructor: compute label combination.
  LabelCombinationLocalShapeBasedAveraging( const UniformVolume::SmartConstPtr targetImage ) : Superclass( targetImage ) {}
  
  /// Add an atlas (pair of reformatted, target-matched intensity image and label map).
  void AddAtlas( const UniformVolume::SmartConstPtr image, const UniformVolume::SmartConstPtr atlas );

  /// Get resulting combined segmentation.
  TypedArray::SmartPtr GetResult() const;  

private:
  /// Compute result for a region.
  void ComputeResultForRegion( const Self::TargetRegionType& region, TypedArray& result ) const;

  /// Signed distance maps for the atlas label maps.
  std::vector<UniformVolume::SmartConstPtr> m_AtlasDMaps;
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationLocalShapeBasedAveraging_h_included_
