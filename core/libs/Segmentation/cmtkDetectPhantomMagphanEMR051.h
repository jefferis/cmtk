/*
//
//  Copyright 2012 SRI International
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

#ifdef CMTK_USE_FFTW
#  include <Segmentation/cmtkSphereDetectionMatchedFilterFFT.h>
#endif

#include <vector>

namespace
cmtk
{

/** \addtogroup Segmentation */
//@{

/** Class for detecting landmark locations of the Magphan EMR051 structural imaging phantom.
 */
class DetectPhantomMagphanEMR051
{
public:
  /// Constructor: compute registrations.
  DetectPhantomMagphanEMR051( UniformVolume::SmartConstPtr& phantomImage );
  
  /// Get landmark coordinates.
  std::vector<UniformVolume::SpaceVectorType> GetLandmarks();

private:
  /// Image of the phantom.
  UniformVolume::SmartConstPtr m_PhantomImage;

  /// The sphere detector.
#ifdef CMTK_USE_FFTW
  SphereDetectionMatchedFilterFFT m_SphereDetector;
#endif

  /// Find a number of spheres of equal size.
  void FindSpheres( std::vector<cmtk::UniformVolume::SpaceVectorType>::iterator dest, const int nSpheres, const Types::Coordinate radius, UniformVolume::SmartPtr& excludeMask );
};

//@}

} // namespace cmtk
