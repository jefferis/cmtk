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

#ifndef __cmtkDeblurringVolumeReconstruction_h_included_
#define __cmtkDeblurringVolumeReconstruction_h_included_

#include <cmtkconfig.h>

#include <cmtkInverseInterpolationVolumeReconstructionBase.h>

#include <cmtkAffineRegistration.h>
#include <cmtkAffineXform.h>
#include <cmtkUniformVolume.h>

#include <ap.h>
#include <vector>

namespace
cmtk
{

/** \addtogroup Recon */
//@{
/** Class for volume reconstruction using joint iterative deblurring.
 *
 *\author Torsten Rohlfing and Michael P. Hasak
 */
template<class TPSF>
class DeblurringVolumeReconstruction :
  /// Inherit from non-templated base class.
  public InverseInterpolationVolumeReconstructionBase
{
public:
  /// This class.
  typedef DeblurringVolumeReconstruction<TPSF> Self;

  /// Superclass.
  typedef InverseInterpolationVolumeReconstructionBase Superclass;

  /** Constructor for interleaved image motion correction.
   * Take original image. Set interleaved image count and stacking axis. Construct separate 3D image stacks for interleaved passes. Allocate corrected image.
   *\param originalImage Smart pointer to the original image with motion artifacts.
   *\param numberOfPasses The number of interleaved passes, i.e., the number of pass images that comprise the final image.
   *\param interleaveAxis Between-slice axis of the interleaved acquisition.
   */
  DeblurringVolumeReconstruction( const UniformVolume* originalImage, const size_t numberOfPasses, const int interleaveAxis, const Types::Coordinate psfScale = 1.0 )
    : InverseInterpolationVolumeReconstructionBase( originalImage, numberOfPasses, interleaveAxis )
  { 
    this->m_FunctionAndGradient = new typename Self::FunctionAndGradient( this );
    
    Vector3D scale( originalImage->m_Delta );
    scale *= psfScale;
    scale.XYZ[interleaveAxis] *= numberOfPasses;
    for ( size_t i = 0; i < numberOfPasses; ++i )
      {
      typename TPSF::SmartPtr psf( new TPSF( scale ) );
      this->m_PassImagePSF.push_back( psf );
      }
  }

  /** Constructor for general volume reconstruction from multiple acquired images.
   */
  DeblurringVolumeReconstruction( const UniformVolume* reconstructionGrid, std::vector<UniformVolume::SmartPtr>& images, const Vector3D& psfSize )
    : InverseInterpolationVolumeReconstructionBase( reconstructionGrid, images )
  { 
    this->m_FunctionAndGradient = new typename Self::FunctionAndGradient( this );

    for ( size_t i = 0; i < images.size(); ++i )
      {
      const Vector3D scale( psfSize );
      typename TPSF::SmartPtr psf( new TPSF( scale ) );
      this->m_PassImagePSF.push_back( psf );
      }
  }
  
private:
  /// Scale factors for acquired images point spread functions.
  std::vector<typename TPSF::SmartPtr> m_PassImagePSF;

  /// Blur corrected image into pass images.
  void Blur( const ap::real_1d_array& reconstructedPixelArray );

  /// Compute gradient of approximation error w.r.t. corrected image pixels.
  void ComputeErrorGradientImage( ap::real_1d_array& g );

  /** Get a bounding box of the transformed pass-image pixel neighborhood.
   * (Transformed from pass-image space to corrected-image space)
   */
  void GetBoundingBoxOfXformedPassNeighborhood
  ( int* region, const UniformVolume* correctedImage, const Vector3D& currentPassVoxel, const TPSF* psf, const AffineXform* passToCorrectedXform, const int* correctedImageDims ) const;

  /// Glue class for function and gradient evaluation.
  class FunctionAndGradient : 
    /// Inherit from virtual base class.
    public ap::FunctionAndGradient
  {
  public:
    /// Function class type.
    typedef DeblurringVolumeReconstruction FunctionType;

    /// Constructor.
    FunctionAndGradient( FunctionType* function )
    {
      this->m_Function = function;
    }

    /// Evaluate function and gradient.
    virtual void Evaluate( const ap::real_1d_array& x, ap::real_value_type& f, ap::real_1d_array& g );

  private:
    /// Pointer to actual function class.
    FunctionType* m_Function;
  };

  /// Give function and gradient evaluator private access.
  friend class Self::FunctionAndGradient;
};

#include <cmtkDeblurringVolumeReconstruction.txx>

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDeblurringVolumeReconstruction_h_included_

