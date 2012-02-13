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

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Smart pointer to const to this class.
  typedef SmartConstPointer<Self> SmartConstPtr;

  /// Constructor: compute label combination.
  LabelCombinationLocalWeighting( const UniformVolume::SmartConstPtr targetImage ) : 
    m_TargetImage( targetImage ),
    m_PatchRadius( UniformVolume::IndexType::Init( 1 ) ),
    m_PatchRadiusPlusOne( UniformVolume::IndexType::Init( 2 ) ),
    m_SearchRegion( UniformVolume::IndexType( UniformVolume::IndexType::Init( 0 ) ),
		    UniformVolume::IndexType( UniformVolume::IndexType::Init( 1 ) ) )
  {}
  
  /// Add an atlas image (reformatted, target-matched intensity image).
  void AddAtlasImage( const UniformVolume::SmartConstPtr image );

  /** Exclude global outliers.
   * Detect atlases with abnormally low correlation between reformatted atlas and target image,
   * then delete these atlases. Outliers are defined as NCC below Q1-1.5*(Q3-Q1), where Q1 is
   * the 25th percentile of NCC between atlas and target over all atlases, Q3 is the 75th 
   * percentile.
   */
  void ExcludeGlobalOutliers();

  /// Set patch radius.
  void SetPatchRadius( const size_t radius )
  {
    this->m_PatchRadius = UniformVolume::IndexType( UniformVolume::IndexType::Init( radius ) );
    this->m_PatchRadiusPlusOne = UniformVolume::IndexType( UniformVolume::IndexType::Init( radius+1 ) );
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

  /// Image patch radius in pixels (x,y,z) plus one added to each dimension.
  UniformVolume::IndexType m_PatchRadiusPlusOne;

  /// Patch search region in pixels (x,y,z).
  UniformVolume::RegionType m_SearchRegion;

  /** Delete atlas with given index. 
   * Derived classes may need to overload this to make sure additional  atlas components (e.g.,
   * distance map, label map) are also properly deleted.
   */
  virtual void DeleteAtlas( const size_t i )
  {
    this->m_AtlasImages.erase( this->m_AtlasImages.begin() + i );
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkLabelCombinationLocalWeighting_h_included_
