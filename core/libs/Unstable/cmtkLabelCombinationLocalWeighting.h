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

#ifndef __cmtkLabelCombinationLocalWeighting_h_included_
#define __cmtkLabelCombinationLocalWeighting_h_included_

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

/** Base class for segmentation combination using local weighting.
 */
class LabelCombinationLocalWeighting
{
public:
  /// This class.
  typedef LabelCombinationLocalWeighting Self;

  /// Constructor: compute label combination.
  LabelCombinationLocalWeighting( const UniformVolume::SmartConstPtr targetImage ) : 
    m_TargetImage( targetImage ),
    m_PatchRadius( UniformVolume::IndexType::Init( 1 ) ),
    m_SearchRegion( UniformVolume::IndexType( UniformVolume::IndexType::Init( 0 ) ),
		    UniformVolume::IndexType( UniformVolume::IndexType::Init( 1 ) ) )
  {}
  
  /// Add an atlas (pair of reformatted, target-matched intensity image and label map).
  void AddAtlasImage( const UniformVolume::SmartConstPtr image );

  /// Set patch radius.
  void SetPatchRadius( const size_t radius )
  {
    this->m_PatchRadius = UniformVolume::IndexType( UniformVolume::IndexType::Init( radius ) );
  }

  /// Set patch radius.
  void SetSearchRadius( const size_t radius )
  {
    this->m_SearchRegion.From() = UniformVolume::IndexType( UniformVolume::IndexType::Init( -radius ) );
    this->m_SearchRegion.To() = UniformVolume::IndexType( UniformVolume::IndexType::Init( radius+1 ) );
  }
  
  /// Get resulting combined segmentation.
  virtual TypedArray::SmartPtr GetResult() const = 0;
  
protected:
  /// Target image region type.
  typedef UniformVolume::RegionType TargetRegionType;

  /// The target image.
  UniformVolume::SmartConstPtr m_TargetImage;
  
  /// Vector of target-matched atlas images.
  std::vector<UniformVolume::SmartConstPtr> m_AtlasImages;

  /// Image patch radius in pixels (x,y,z).
  UniformVolume::IndexType m_PatchRadius;

  /// Patch search region in pixels (x,y,z).
  UniformVolume::RegionType m_SearchRegion;
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationLocalWeighting_h_included_
