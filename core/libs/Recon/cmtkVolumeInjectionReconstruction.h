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

#ifndef __cmtkVolumeInjectionReconstruction_h_included_
#define __cmtkVolumeInjectionReconstruction_h_included_

#include <cmtkconfig.h>

#include <cmtkAffineRegistration.h>
#include <cmtkAffineXform.h>
#include <cmtkXform.h>
#include <cmtkUniformVolume.h>

#include <vector>
#include <ap.h>

namespace
cmtk
{

/** \addtogroup Recon */
//@{

/** Class for volume reconstruction using volume injection.
 *
 * This class implements some side effects that are not strictly part of the volume
 * injection algorithm, but which supply functionality to the derived
 * igsInverseInterpolationVolumeReconstructionBase class. These include computation
 * of regional minimum/maximum intensity ranges, between-pass registration, and
 * KLD metric computation.
 *
 *\author Torsten Rohlfing
 */
class VolumeInjectionReconstruction
{
public:
  /// This class.
  typedef VolumeInjectionReconstruction Self;

  /** Constructor from single interleaved image.
   * Take original image. Set interleaved image count and stacking axis. Construct separate 3D image stacks for interleaved passes. Allocate corrected image.
   *\param originalImage Smart pointer to the original image with motion artifacts.
   *\param numberOfPasses The number of interleaved passes, i.e., the number of pass images that comprise the final image.
   *\param interleaveAxis Between-slice axis of the interleaved acquisition.
   */
  VolumeInjectionReconstruction( const UniformVolume* originalImage, const int numberOfPasses, const int interleaveAxis );

  /** Constructor for general volume reconstruction from multiple acquired images.
   */
  VolumeInjectionReconstruction( const UniformVolume* reconstructionGrid, std::vector<UniformVolume::SmartPtr>& images );

  /// Virtual destructor stub.
  virtual ~VolumeInjectionReconstruction() {}

  /** Static helper function: guess interleaved axis.
   * Basically, we assume that images are acquired as interleaved stacks of
   * square 2D images with pixel size different from inter-slice spacing.
   * So we guess that the interleaved axis is the one that does not match the
   * other two axes' image dimensions, and if all three are the same, the
   * one that doesn't match their spacing (delta).
   */
  static int GuessInterleaveAxis( 
    const UniformVolume* image, //!< The interleaved image.
    const int defaultAxis = 2 //!< In case all guessing fails, this is the default axis we return.
    );

  /** Compute transformations between the reference image grid and the original pass images.
   * The resulting transformations are stored in the m_TransformationsToPassImages vector.
   * If a high-resolution reference image is set in m_ReferenceImage, then all subimages are registered to it.
   * Otherwise, the 0th subimage is used as the reference and the remaining subimages are registered to it.
   * In this case, the transformation for the 0th subimage is set as the identity transformation.
   *\param registrationMetric Similarity metric for registration of the passes to the reference image.
   */
  void ComputeTransformationsToPassImages( const int registrationMetric = 0 );

  /// Set transformations to pass images externally (e.g., imported from disk).
  void SetTransformationsToPassImages( std::vector<Xform::SmartPtr>& transformations )
  {
    this->m_TransformationsToPassImages = transformations;
  }

  /// Get transformation to one pass image.
  Xform::SmartPtr& GetTransformationToPassImage( const size_t passIdx )
  {
    if ( passIdx < this->m_TransformationsToPassImages.size() )
      return this->m_TransformationsToPassImages[passIdx];
    else
      return Xform::SmartPtr::Null;
  }

  /// Get transformation to one pass image.
  std::vector<Xform::SmartPtr>& GetTransformationsToPassImages()
  {
    return this->m_TransformationsToPassImages;
  }

  /// Create initial approximation using isotropic volume injection.
  void VolumeInjectionIsotropic
  ( const Types::Coordinate kernelSigma, //!< Gaussian kernel sigma (standard deviation) parameter
    const Types::Coordinate kernelRadius ); //!< Gaussian kernel cut-off radius.

  /// Create initial approximation using anisotropic volume injection.
  void VolumeInjectionAnisotropic
  ( const Types::Coordinate kernelSigmaFactor, //!< Gaussian kernel sigma (standard deviation) factor (multiple of per-dimension pass image spacing)
    const Types::Coordinate kernelRadiusFactor ); //!< Gaussian kernel cut-off radius factor (multiple of per-dimension pass image spacing)

  /// Returns the corrected image.
  UniformVolume::SmartPtr& GetCorrectedImage();
  
  /// Set optional separate reference image for motion parameter estimation.
  void SetReferenceImage( UniformVolume::SmartPtr& referenceImage );  

  /** Set pass weight.
   * Each pass weight should be between 0 and 1. If the weight for a pass is zero, then that
   * pass is effectively excluded from the reconstruction. This can be useful if one of the
   * passes shows severe within-pass motion artifacts that would otherwise disturb the
   * across-pass correction.
   *
   * By default, all pass weights are set to 1, i.e., all passes contribute equally.
   */
  void SetPassWeight( const size_t pass, const Types::Coordinate weight )
  {
    this->m_PassWeights[pass] = weight;
  }

  /// Get Kullback-Leibler Divergence between intensity distributions in original and corrected image.
  ap::real_value_type GetOriginalToCorrectedImageKLD( const ap::real_1d_array& x );
  
protected:
  /// Number of interleaved passes.
  int m_NumberOfPasses;

  /// Relative weights of the passes in the correction; can be used to underweight or even exclude passes.
  std::vector<Types::Coordinate> m_PassWeights;

  /// Original volume pixel intensity minimum.
  Types::DataItem m_OriginalImageMax;

  /// Original volume pixel intensity minimum.
  Types::DataItem m_OriginalImageMin;

  /// Original pass images.
  std::vector<UniformVolume::SmartPtr> m_OriginalPassImages;

  /// Histogram type.
  typedef Histogram<double> HistogramType;

  /// Original image histogram.
  HistogramType::SmartPtr m_OriginalImageHistogram;

  /// Corrected image histogram.
  HistogramType::SmartPtr m_CorrectedImageHistogram;

  /// Original image intensity noise kernel.
  std::vector<HistogramType::BinType> m_OriginalImageIntensityNoiseKernel;

  /// Optional high-resolution non-interleaved reference image.
  UniformVolume::SmartPtr m_ReferenceImage;

  /// Affine transformations that map FROM the corrected image TO each of the subimages.
  std::vector<Xform::SmartPtr> m_TransformationsToPassImages;

  /// Developing corrected image.
  UniformVolume::SmartPtr m_CorrectedImage;

  /// Corrected image Laplacian.
  std::vector<ap::real_value_type> m_CorrectedImageLaplacians;  

  /** Compute norm of the corrected image Laplacian.
   * Side effect: this function first computes the Laplacian image, which is stored in
   * m_CorrectedImageLaplacian for use in the AddLaplacianGradientImage function.
   */
  ap::real_value_type ComputeCorrectedImageLaplacianNorm( 
    const ap::real_1d_array& correctedImagePixels //!< Current vector of corrected image pixels.
    );
  
  /// Add weighted gradient image of Laplacian to already computed cost function gradient.
  void AddLaplacianGradientImage( 
    ap::real_1d_array& g, 
    const ap::real_1d_array& correctedImagePixels,
    const ap::real_value_type weight ) const;

  /// Maximum neighborhood pixel values in the corrected image.
  ap::real_1d_array m_NeighorhoodMaxPixelValues;
  
  /// Minimum neighborhood pixel values in the corrected image.
  ap::real_1d_array m_NeighorhoodMinPixelValues;
  
private:
  /// Number of histogram bins for image entropy estimation.
  static const unsigned int NumberOfHistogramBins = 64;

  /// Setup kernels and histograms for image entropy estimation.
  void SetupHistogramKernels( const TypedArray* originalData );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVolumeInjectionReconstruction_h_included_

