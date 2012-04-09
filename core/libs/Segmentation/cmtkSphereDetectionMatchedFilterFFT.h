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

#include <fftw3.h>

namespace
cmtk
{

/** Detect spheres in an image using FFT-based matched filter.
 *\note This class requires CMTK to be configured with FFTW3 support ("CMTK_USE_FFTW" CMake option).
 */
class SphereDetectionMatchedFilterFFT
{
public:
  /// Constructor: initialize FFTW plans and compute image FT.
  SphereDetectionMatchedFilterFFT( const UniformVolume& image );

  /// Destructor: destroy FFTW plans and de-allocate transformed image memories.
  virtual ~SphereDetectionMatchedFilterFFT();

  /// Get image filtered with spherical matched filter kernel.
  cmtk::TypedArray::SmartPtr GetFilteredImageData( const Types::Coordinate sphereRadius /*!< Radius of detected spheres in world coordinate units (e.g., mm) */, 
						   const int marginWidth = 1 /*!< Half width of the filter margin in pixels: positive filter coefficients in a band of this width inside radius, negative coeffiecients outside radius.*/ );

private:
  /// Image number of pixels.
  size_t m_NumberOfPixels;

  /// Image dimensions.
  DataGrid::IndexType m_ImageDims;

  /// Image pixel size.
  UniformVolume::SpaceVectorType m_PixelSize;

  /// The Fourier-transformed image.
  fftw_complex* m_ImageFT;

  /// The Fourier-transformed matched filter.
  fftw_complex* m_FilterFT;

  /// The filter FFT plan.
  fftw_plan m_PlanFilter;

  /// The backward (filtered data to space domain) FFT plan.
  fftw_plan m_PlanBackward;

  /** Make the filter kernel.
   *\return Number of non-zero elements in the filter kernel.
   */
  size_t MakeFilter( const Types::Coordinate sphereRadius, const int marginWidth );
};

} // namespace cmtk
