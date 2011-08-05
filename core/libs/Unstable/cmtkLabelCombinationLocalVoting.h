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

#ifndef __cmtkLabelCombinationLocalVoting_h_included_
#define __cmtkLabelCombinationLocalVoting_h_included_

#include <cmtkconfig.h>

#include <System/cmtkSmartPtr.h>
#include <System/cmtkSmartConstPtr.h>

#include <Base/cmtkUniformVolume.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Segmentation combination by locally-weighted label voting.
 *\attention All labels must be between 0 and 255.
 */
class LabelCombinationLocalVoting
{
public:
  /// Constructor: compute label combination.
  LabelCombinationLocalVoting( const UniformVolume::SmartConstPtr targetImage ) : m_TargetImage( targetImage ) {}
  
  /// Add an atlas (pair of reformatted, target-matched intensity image and label map).
  void AddAtlas( const UniformVolume::SmartConstPtr image, const UniformVolume::SmartConstPtr atlas );

  /// Get resulting combined segmentation.
  TypedArray::SmartPtr GetResult() const;
  
private:
  /// The target image.
  UniformVolume::SmartConstPtr m_TargetImage;
  
  /// Vector of target-matched atlas images.
  std::vector<UniformVolume::SmartConstPtr> m_AtlasImages;

  /// Vector of target-matched atlas label maps.
  std::vector<UniformVolume::SmartConstPtr> m_AtlasLabels;

  /// Neighborhood radius in pixels (x,y,z).
  UniformVolume::IndexType m_NeighborhoodRadius;
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationLocalVoting_h_included_