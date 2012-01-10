/*
//
//  Copyright 2011, 2012 SRI International
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

#ifndef __cmtkEchoPlanarUnwarpFunctional_h_included_
#define __cmtkEchoPlanarUnwarpFunctional_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>

#include <Numerics/ap.h>
#include <Numerics/lbfgsb.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for computing a const function for unwarping Echo Planar Images.
 *\see D. Holland, J. M. Kuperman, and A. M. Dale, "Efficient correction of inhomogeneous static magnetic field-induced distortion in Echo Planar Imaging," NeuroImage, vol. 50, no. 1, pp. 175-183, 2010.
 */
class EchoPlanarUnwarpFunctional
{
public:
  /// This class.
  typedef EchoPlanarUnwarpFunctional Self;

  /// Sinc image interpolation kernel radius.
  static const int InterpolationKernelRadius; 

  /// Constructor.
  EchoPlanarUnwarpFunctional( UniformVolume::SmartConstPtr& imageFwd /*!< "Forward direction" EPI */, 
			      UniformVolume::SmartConstPtr& imageRev /*!< "Reverse" direction EPI */, 
			      const byte phaseEncodeDirection /*!< Phase encoding direction (image coordinate axis) */ );

  /** Set image smoothing kernel width and create smoothed images.
   *\note Smoothing is applied only in phase encoding direction because the deformation field is restricted to this direction, and the optimization can thus not
   * get "caught" in local minima in the other two directions.
   */
  void SetSmoothingKernelWidth ( const Units::GaussianSigma& sigma, /*!< Kernel parameter "sigma" (standard deviation). Setting kernel width to "0" simply uses original images. */
				 const Types::Coordinate maxError = 1e-5 /*!< Maximum approximation error: the kernel is truncated when it falls below this threshold */ );

  /// Return either first or second corrected image.
  UniformVolume::SmartPtr GetCorrectedImage( const byte idx = 0 )
  {
    UniformVolume::SmartPtr correctedImage( this->m_ImageFwd->CloneGrid() );

    const std::vector<Types::DataItem>& srcImage = ( idx == 0 ) ? this->m_CorrectedImageFwd : this->m_CorrectedImageRev;

    correctedImage->CreateDataArray( TYPE_FLOAT );
    for ( size_t px = 0; px < this->m_ImageFwd->GetNumberOfPixels(); ++px )
      {
      correctedImage->SetDataAt( srcImage[px], px );
      }

    return correctedImage;
  }
  
  /// Optimize unwarping deformation using L-BFGS optimizer.
  void Optimize( const int numberOfIterations, const Units::GaussianSigma& smoothMax, const Units::GaussianSigma& smoothMin, const Units::GaussianSigma& smoothDiff );
  
  /** Set smoothness constraint weight.
   * This is \lambda_2 in Eq. (9) in Holland et al. NIMG 2010.
   */
  void SetSmoothnessConstraintWeight( const Types::Coordinate weight )
  {
    this->m_SmoothnessConstraintWeight = weight;
  }
  
  /** Set folding constraint weight.
   * This is based on the log() of the local 1D Jacobians of the deformation field.
   */
  void SetFoldingConstraintWeight( const Types::Coordinate weight )
  {
    this->m_FoldingConstraintWeight = weight;
  }
  
private:
  /// Unwarped image grid.
  UniformVolume::SmartPtr m_ImageGrid;

  /// "Forward" phase encoding image.
  UniformVolume::SmartConstPtr m_ImageFwd;

  /// "Reverse" phase encoding image.
  UniformVolume::SmartConstPtr m_ImageRev;

  /// Smoothed "forward" phase encoding image.
  UniformVolume::SmartConstPtr m_SmoothImageFwd;

  /// Smoothed "reverse" phase encoding image.
  UniformVolume::SmartConstPtr m_SmoothImageRev;

  /// Phase encoding direction.
  byte m_PhaseEncodeDirection;

  /** Readout direction.
   * This is a bit misnamed - really this is the coordinate dimension, excluding the phase encoding direction,
   * with the maximum number of pixels. The purpose of this is for a maximally efficient partition for SMP
   * computation. So while the name may not always fit the determined value, and there may be cases when this is
   * selected as the slice direction, it does serve its purpose.
   */
  byte m_ReadoutDirection;

  /// Smoothness constraint weight.
  Types::Coordinate m_SmoothnessConstraintWeight;

  /// Folding constraint weight.
  Types::Coordinate m_FoldingConstraintWeight;

  /// 1D deformation map along phase encoding direction.
  ap::real_1d_array m_Deformation;

  /// 1D intensity gradient of "forward" image.
  std::vector<Types::DataItem> m_GradientImageFwd;

  /// 1D intensity gradient of "reverse" image.
  std::vector<Types::DataItem> m_GradientImageRev;
  
  /// Deformed "forward" image.
  std::vector<Types::DataItem> m_UnwarpImageFwd;

  /// Deformed "reverse" image.
  std::vector<Types::DataItem> m_UnwarpImageRev;

  /// Deformed and intensity-corrected "forward" image.
  std::vector<Types::DataItem> m_CorrectedImageFwd;

  /// Deformed and intensity-corrected "reverse" image.
  std::vector<Types::DataItem> m_CorrectedImageRev;

  /// Compute 1D intensity gradient image.
  void MakeGradientImage( const ap::real_1d_array& x /*!< Current deformation parameter vector.*/, 
			  const int direction /*!< +1 = forward image, -1 = reverse image */,
			  const UniformVolume& sourceImage /*!< Image to compute gradient for.*/, 
			  std::vector<Types::DataItem>& gradientImageData /*!< Reference to data array for computed gradient data.*/ );

  /// Compute deformed image.
  void ComputeDeformedImage( const ap::real_1d_array& x /*!< Current deformation parameter vector.*/, 
			     int direction /*!< Deformation direction - 1 computes unwarped "forward" image, -1 computed unwarped "reverse" image.*/,
			     const UniformVolume& sourceImage /*!< Undeformed input image.*/, 
			     std::vector<Types::DataItem>& targetUnwarpData /*!< Reference to deformed output image data.*/,
			     std::vector<Types::DataItem>& targetCorrectedData /*!< Reference to deformed and intensity-corrected output image data.*/ );

  /// 1D sinc interpolation
  Types::DataItem Interpolate1D( const UniformVolume& sourceImage /*!< Image to interpolate from */, 
				 const FixedVector<3,int>& baseIdx /*!< Grid base index - this is the grid cell where the interpolation kernel is anchored. */, 
				 const Types::Coordinate relative /*!< Relative position of interpolation location in grid cell, relative to phase encoding direction. */ ) const;

  /** Get partial image deformation Jacobian.
   *\return The forward image Jacobian can be computed from this as Jfwd = 1+Jpartial, the reverse Jacobian as Jrev = 1-Jpartial.
   */
  Types::Coordinate GetPartialJacobian( const ap::real_1d_array& x /*!< Current deformation parameter vector.*/, 
					const FixedVector<3,int>& baseIdx /*!< Grid base index - this is the grid cell where the differential operator stencil is anchored. */ ) const;

  /// Glue class for function and gradient evaluation.
  class FunctionAndGradient : 
    /// Inherit from virtual base class.
    public ap::FunctionAndGradient
  {
  public:
    /// Function class type.
    typedef EchoPlanarUnwarpFunctional FunctionType;

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

//@}

} // namespace cmtk

#endif // #ifndef __cmtkEchoPlanarUnwarpFunctional_h_included_
