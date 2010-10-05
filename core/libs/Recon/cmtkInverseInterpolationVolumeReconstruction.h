/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkInverseInterpolationVolumeReconstruction_h_included_
#define __cmtkInverseInterpolationVolumeReconstruction_h_included_

#include <cmtkconfig.h>

#include <Recon/cmtkInverseInterpolationVolumeReconstructionBase.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkUniformVolume.h>

#include <Registration/cmtkAffineRegistration.h>

#include <System/cmtkProgress.h>

#include "Numerics/ap.h"

#include <vector>

namespace
cmtk
{

/** \addtogroup Recon */
//@{

/** Class for volume reconstruction using inverse interpolation.
 *\author Torsten Rohlfing
 */
template<class TInterpolator> 
class InverseInterpolationVolumeReconstruction :
  /// Inherit from non-templated base class.
  public InverseInterpolationVolumeReconstructionBase
{
public:
  /// This class.
  typedef InverseInterpolationVolumeReconstruction<TInterpolator> Self;

  /// Superclass.
  typedef InverseInterpolationVolumeReconstructionBase Superclass;

  /** Constructor for interleaved image motion correction.
   * Take original image. Set interleaved image count and stacking axis. Construct separate 3D image stacks for interleaved passes. Allocate corrected image.
   *\param originalImage Smart pointer to the original image with motion artifacts.
   *\param numberOfPasses The number of interleaved passes, i.e., the number of pass images that comprise the final image.
   *\param interleaveAxis Between-slice axis of the interleaved acquisition.
   */
  InverseInterpolationVolumeReconstruction( const UniformVolume* originalImage, const int numberOfPasses, const int interleaveAxis )
    : InverseInterpolationVolumeReconstructionBase( originalImage, numberOfPasses, interleaveAxis )
  { 
    this->m_FunctionAndGradient = new typename Self::FunctionAndGradient( this );
  }

  /** Constructor for general volume reconstruction from multiple acquired images.
   */
  InverseInterpolationVolumeReconstruction( const UniformVolume* reconstructionGrid, std::vector<UniformVolume::SmartPtr>& images )
    : InverseInterpolationVolumeReconstructionBase( reconstructionGrid, images )
  { 
    this->m_FunctionAndGradient = new typename Self::FunctionAndGradient( this );
  }
  
  /// Destructor: delete function evaluator object.
  virtual ~InverseInterpolationVolumeReconstruction()
  {
    delete this->m_FunctionAndGradient;
  }

private:
  /// Interpolates subimages from the existing corrected image.
  void Interpolation( const ap::real_1d_array& reconstructedPixelArray );
  
  /// Compute gradient of approximation error w.r.t. corrected image pixels.
  void ComputeErrorGradientImage( ap::real_1d_array& g );

  /// Get pass image region that contains all pixels dependent on currently considered corrected image pixel.
  void GetPassImageDependentPixelRegion
  ( int* region, const UniformVolume* correctedImage, const int* currentCorrectedGridPoint, 
    const UniformVolume* passImage, const AffineXform* transformationToPassImage, const DataGrid::IndexType& passImageDims );

  /// Glue class for function and gradient evaluation.
  class FunctionAndGradient : 
    /// Inherit from virtual base class.
    public ap::FunctionAndGradient
  {
  public:
    /// Function class type.
    typedef InverseInterpolationVolumeReconstruction<TInterpolator> FunctionType;

    /// Constructor.
    FunctionAndGradient( FunctionType* function )
    {
      this->m_Function = function;
    }

    /// Evaluate function and gradient.
    virtual void Evaluate( const ap::real_1d_array& x, ap::real_value_type& f, ap::real_1d_array& g );

    /// Get notified when L-BFGS-B goes into next iteration.
    virtual void NextIteration( const int iteration )
    {
      Progress::SetProgress( iteration );
    }

  private:
    /// Pointer to actual function class.
    FunctionType* m_Function;
  };

  /// Give function and gradient evaluator private access.
  friend class Self::FunctionAndGradient;
};

//@}

} // namespace cmtk

#include "cmtkInverseInterpolationVolumeReconstruction.txx"

#endif // #ifndef __cmtkInverseInterpolationVolumeReconstruction_h_included_

