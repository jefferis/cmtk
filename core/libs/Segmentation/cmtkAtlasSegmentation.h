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

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkSplineWarpXform.h>

namespace
cmtk
{

/** Class for integrated atlas-based segmentation.
 * This class encapsulates the following stages of atlas-based image segmentation: 1) linear image-to-atlas registration, 2) non-linear
 * image-to-atlas registration, 3) label map reformatting.
 */
class AtlasSegmentation
{
public:
  /// Constructor: compute registrations.
  AtlasSegmentation( UniformVolume::SmartPtr& targetImage, UniformVolume::SmartPtr& atlasImage, UniformVolume::SmartPtr& atlasLabels );

  /// Get affine transformation.
  AffineXform::SmartPtr& GetAffineXform()
  {
    if ( ! this->m_AffineXform )
      this->RegisterAffine();
    return this->m_AffineXform;
  }
  
  /// Get nonrigid transformation.
  WarpXform::SmartPtr GetWarpXform()
  {
    if ( ! this->m_WarpXform )
      this->RegisterSpline();
    return this->m_WarpXform;
  }
  
  /// Get nonrigid spline transformation.
  SplineWarpXform::SmartPtr GetSplineWarpXform()
  {
    return SplineWarpXform::SmartPtr::DynamicCastFrom( this->GetWarpXform() );
  }
  
  /// Get reformatted label map.
  UniformVolume::SmartPtr& GetLabelMap()
  {
    if ( ! this->m_LabelMap )
      this->ReformatLabels();
    return this->m_LabelMap;
  }

  /// Set fast flag.
  void SetFast( const bool fast )
  {
    this->m_Fast = fast;
  }

private:
  /// Flag for "fast" computation.
  bool m_Fast;

  /// Target image.
  UniformVolume::SmartPtr m_TargetImage;

  /// Atlas image.
  UniformVolume::SmartPtr m_AtlasImage;

  /// Atlas labels.
  UniformVolume::SmartPtr m_AtlasLabels;

  /// Affine registration transformation.
  AffineXform::SmartPtr m_AffineXform;
  
  /// Compute affine registration.
  void RegisterAffine();

  /// Nonrigid, B-spline transformation.
  WarpXform::SmartPtr m_WarpXform;

  /// Compute spline registration.
  void RegisterSpline();

  /// Output label map.
  UniformVolume::SmartPtr m_LabelMap;

  /// Compute label map.
  void ReformatLabels();
};

} // namespace cmtk
