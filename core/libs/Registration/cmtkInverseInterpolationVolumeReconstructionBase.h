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

#ifndef __cmtkInverseInterpolationVolumeReconstructionBase_h_included_
#define __cmtkInverseInterpolationVolumeReconstructionBase_h_included_

#include <cmtkconfig.h>

#include <cmtkVolumeInjectionReconstruction.h>

#include <cmtkAffineRegistration.h>
#include <cmtkAffineXform.h>
#include <cmtkUniformVolume.h>

#include <vector>

#include <ap.h>
#include <lbfgsb.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for volume reconstruction using inverse interpolation.
 *
 * This class is specifically designed to use the L-BFFS-B Optimizer from alglib.net.
 *
 *\see http://www.alglib.net/optimization/lbfgsb.php
 *
 *\author Torsten Rohlfing
 */
class InverseInterpolationVolumeReconstructionBase :
  /// Inherit from volume injection class.
  public VolumeInjectionReconstruction
{
public:
  /// This class.
  typedef InverseInterpolationVolumeReconstructionBase Self;

  /// Parent class.
  typedef VolumeInjectionReconstruction Superclass;

  /** Constructor for interleaved image motion correction.
   * Take original image. Set interleaved image count and stacking axis. Construct separate 3D image stacks for interleaved passes. Allocate corrected image.
   *\param originalImage Smart pointer to the original image with motion artifacts.
   *\param numberOfPasses The number of interleaved passes, i.e., the number of pass images that comprise the final image.
   *\param interleaveAxis Between-slice axis of the interleaved acquisition.
   */
  InverseInterpolationVolumeReconstructionBase( const UniformVolume* originalImage, const int numberOfPasses, const int interleaveAxis );

  /** Constructor for general volume reconstruction from multiple acquired images.
   */
  InverseInterpolationVolumeReconstructionBase( const UniformVolume* reconstructionGrid, std::vector<UniformVolume::SmartPtr>& images );

  /// Virtual destructor stub.
  virtual ~InverseInterpolationVolumeReconstructionBase() {}

  /** Compute approximation error.
   * This function computed the difference images between the original and the interpolated pass images,
   * and it also computes from these the mean squared and the maximum approximation errors.
   *
   *\return The approximation error, i.e., the sum of squared differences of the intensities
   * in the original pass images and the corresponding intensities in the images interpolated
   * from the current estimated corrected image.
   */
  double ComputeApproximationError();
  
  /// Returns the reconstructed image wioth the lowest maximum error.
  UniformVolume::SmartPtr& GetLowestMaxErrorImage()
  {
    return this->m_LowestMaxErrorImage;
  }
  
  /// Set flag for regional intensity trunctation.
  void SetUseRegionalIntensityTruncation( const bool flag = true )
  {
    this->m_RegionalIntensityTruncation = flag;
  }

  /// Set fourth order error flag.
  void SetUseFourthOrderError( const bool flag = true )
  {
    this->m_FourthOrderError = flag;
  }

  /// Set L-norm constraint weight.
  void SetConstraintWeightLNorm( const double weight )
  {
    this->m_ConstraintWeightLNorm = weight;
  }

  /// Optimize approximation error.
  void Optimize( const int numberOfIterations );

  /// Return MSD error.
  double GetMeanSquaredError() const
  {
    return this->m_MeanSquaredError;
  }

  /// Return maximum error.
  double GetMaximumError() const
  {
    return this->m_MaximumError;
  }

protected:
  /// Flag for regional pixel intensity truncation.
  bool m_RegionalIntensityTruncation;

  /// Corrected image with lowest maximum error.
  UniformVolume::SmartPtr m_LowestMaxErrorImage;

  /// Current lowest maximum error.
  double m_LowestMaxError;
  
  /// Interpolated (from current estimate of correcting image) pass images.
  std::vector<UniformVolume::SmartPtr> m_InterpolatedPassImages;
  
  /// Difference Subimages between original and interpolated Subimages, used for error- and gradientcalculation.
  std::vector<UniformVolume::SmartPtr> m_DifferencePassImages;

  /** Flag for fourth order error optimization.
   * If this flag is set, the error optimized is the fourth order error, otherwise the second order (square) error.
   */
  bool m_FourthOrderError;

  /** Constraint weight for Laplacian norm regularizer.
   * Regularization is turned off if weight is less than, or equal to, zero.
   */
  double m_ConstraintWeightLNorm;

  /// Mean squared error of the approximation, i.e., MSD between interpolated and actual pass images.
  double m_MeanSquaredError;

  /// Maximum error of the approximation, i.e., maximum difference between interpolated and actual pass images.
  double m_MaximumError;

  /// Function and gradient evaluator object.
  ap::FunctionAndGradient* m_FunctionAndGradient;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkInverseInterpolationVolumeReconstructionBase_h_included_

