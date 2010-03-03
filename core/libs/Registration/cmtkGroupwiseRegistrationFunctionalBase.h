/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkGroupwiseRegistrationFunctionalBase_h_included_
#define __cmtkGroupwiseRegistrationFunctionalBase_h_included_

#include <cmtkconfig.h>

#include <cmtkFunctional.h>

#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

#include <cmtkUniformVolume.h>
#include <cmtkXform.h>

#include <vector>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Base class for groupwise registration functionals.
 */
class GroupwiseRegistrationFunctionalBase : 
  /** Inherit from abstract functional base class. */
  public Functional
{
public:
  /// Type of parent class.
  typedef Functional Superclass;

  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalBase Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  GroupwiseRegistrationFunctionalBase();

  /// Destructor.
  virtual ~GroupwiseRegistrationFunctionalBase();

  /** Set flag for freeing and rereading images.
   * If registration uses smoothed images, the original data can be freed after smoothing
   * and reread from the file system if needed again. This saves roughly 1/2 of memory
   * allocation.
   */
  virtual void SetFreeAndRereadImages( const bool flag = true )
  {
    this->m_FreeAndRereadImages = flag;
  }

  /// Set flag for repeated histogram-based intensity matching.
  virtual void SetRepeatIntensityHistogramMatching( const bool flag = true )
  {
    this->m_RepeatIntensityHistogramMatching = flag;
  }

  /** Create template grid based on target images.
   * The template image is constructed so that is has the maximum number of
   * pixels of all targets in each dimension, and the minimum isotropic
   * pixel size.
   *
   *\param downsample Downsampling factor. The voxel size in the template
   * image is increased by this factor.
   */
  virtual void CreateTemplateGridFromTargets( const std::vector<UniformVolume::SmartPtr>& targets, const int downsample = 0 );

  /** Create template grid based on geometry.
   */
  virtual void CreateTemplateGrid( const int (&dims)[3], const Types::Coordinate (&deltas)[3] );

  /** Set template grid.
   *\param The template grid that defines size and resolution for the
   *  implicit registration template.
   */
  virtual void SetTemplateGrid( UniformVolume::SmartPtr& templateGrid, const int downsample = 1, const bool useTemplateData = false );

  /** Retrieve the template grid.
   */
  virtual UniformVolume::SmartPtr& GetTemplateGrid() 
  { 
    return this->m_TemplateGrid; 
  }

  /** Retrieve the template grid.
   */
  virtual const UniformVolume* GetTemplateGrid() const
  {
    return this->m_TemplateGrid; 
  }

  /** Set target images.
   *\param tImages Vector of all images to be registered.
   */
  virtual void SetTargetImages( std::vector<UniformVolume::SmartPtr>& tImages );

  /** Get the original target images.
   */
  virtual void GetOriginalTargetImages( std::vector<UniformVolume::SmartPtr>& tImages )
  {
    tImages = this->m_OriginalImageVector;
  }  

  /** Get the original target images.
   */
  virtual std::vector<UniformVolume::SmartPtr>& GetOriginalTargetImages()
  {
    return this->m_OriginalImageVector;
  }  

  /** Get number of target images.
   */
  virtual size_t GetNumberOfTargetImages() const
  {
    return this->m_OriginalImageVector.size();
  }  

  /** Get a smart pointer to one original target image.
   */
  virtual UniformVolume::SmartPtr GetOriginalTargetImage( const size_t imageIdx )
  {
    return this->m_OriginalImageVector[imageIdx];
  }  

  /** Get a constant pointer to one original target image.
   */
  virtual const UniformVolume* GetOriginalTargetImage( const size_t imageIdx ) const
  {
    return this->m_OriginalImageVector[imageIdx];
  }  

  /** Set Gaussian smoothing kernel width for target images.
   * Non-positive values turn off smoothing.
   */
  virtual void SetGaussianSmoothImagesSigma( const Types::Coordinate gaussianSmoothImagesSigma )
  {
    this->m_GaussianSmoothImagesSigma = gaussianSmoothImagesSigma;
  }

  /** Set user-defined image background value to be substituted outside the field of view. */
  virtual void SetUserBackgroundValue( const byte value = 0 )
  {
    this->m_UserBackgroundValue = value;
    this->m_UserBackgroundFlag = true;
  }

  /** Unset user-defined image background value. */
  virtual void UnsetUserBackgroundValue()
  {
    this->m_UserBackgroundFlag = false;
  }

  /// Set flag for zero-sum updates.
  virtual void SetForceZeroSum( const bool forceZeroSum = true )
  {
    this->m_ForceZeroSum = forceZeroSum;
  }

  /// Set count for restricted zero-sum updates.
  virtual void SetForceZeroSumFirstN( const size_t forceZeroSumFirstN )
  {
    this->m_ForceZeroSumFirstN = forceZeroSumFirstN;
    this->m_ForceZeroSum = this->m_ForceZeroSum || (this->m_ForceZeroSumFirstN>0);
  }

  /** Set range of currently active images.
   * The "to" parameter is the index of the last active image plus one, so
   * it can be used directly as the upper bound in a "for" loop.
   */
  virtual void SetActiveImagesFromTo( const size_t from, const size_t to )
  {
    this->m_ActiveImagesFrom = from;
    this->m_ActiveImagesTo = to;
  }

  /** Set range of currently active transformations.
   * The "to" parameter is the index of the last active transformation plus one, so
   * it can be used directly as the upper bound in a "for" loop.
   */
  virtual void SetActiveXformsFromTo( const size_t from, const size_t to )
  {
    this->m_ActiveXformsFrom = from;
    this->m_ActiveXformsTo = to;
  }

  /** Set probabilistic sampling density.
   */
  virtual void SetProbabilisticSampleDensity( const float density )
  {
    this->m_ProbabilisticSampleDensity = density;
  }

  /** Set number of iterations after which probabilistic samples are updated.
   */
  virtual void SetProbabilisticSampleUpdatesAfter( const int iterations )
  {
    this->m_ProbabilisticSampleUpdatesAfter = iterations;
    this->m_ProbabilisticSampleUpdatesSince = 0;
  }

  /** Set transformations.
   */
  void SetXforms( const std::vector<Xform::SmartPtr>& xformVector )
  {
    this->m_XformVector.resize( xformVector.size() );
    for ( size_t i = 0; i < this->m_XformVector.size(); ++i )
      {
      this->m_XformVector[i] = xformVector[i];
      }
  }

  /** Get coordinate transformation for one image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual const Xform* GetGenericXformByIndex( const size_t idx ) const
  {
    return this->m_XformVector[idx];
  }

  /** Get coordinate transformation for one image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual Xform::SmartPtr GetGenericXformByIndex( const size_t idx )
  {
    return this->m_XformVector[idx];
  }

  /** Get coordinate transformation for one active image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual const Xform* GetGenericActiveXformByIndex( const size_t idx ) const
  {
    return this->m_XformVector[idx + this->m_ActiveXformsFrom];
  }

  /** Get coordinate transformation for one active image in the group.
   *\param idx Index of the volume/transformation.
   *\return Transformation for the selected volume.
   */
  virtual Xform::SmartPtr GetGenericActiveXformByIndex( const size_t idx )
  {
    return this->m_XformVector[idx + this->m_ActiveXformsFrom];
  }

  /** Get parameter stepping in milimeters.
   *\param idx Parameter index.
   *\return Step of given parameter that corresponds to 1 mm effective motion.
   */
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    const size_t xformIdx = idx / this->m_ParametersPerXform;
    if ( (xformIdx >= this->m_ActiveXformsFrom) && (xformIdx < this->m_ActiveXformsTo) )
      {
      return this->m_XformVector[xformIdx]->GetParamStep( idx % this->m_ParametersPerXform, this->m_ImageVector[xformIdx]->Size, mmStep );
      }
    else
      {
      return 0.0;
      }
  }
  
  /** Get parameter stepping in milimeters.
   *\param idx Parameter index.
   *\return Step of given parameter that corresponds to 1 mm effective motion.
   */

  /** Return the functional's parameter vector dimension.
   * We assume that all transformations have the same number of parameters.
   * This is true for affine transformations.
   */
  virtual size_t ParamVectorDim() const 
  { 
    return this->m_ParametersPerXform * this->m_XformVector.size(); 
  }
  
  /** Return the number of variable parameters of the transformation.
   *\return This function returns the same value as ParamVectorDim(). 
   *  Non-varying parameters (e.g., rotation centers) are handled via
   *  parameter step values.
   */
  virtual size_t VariableParamVectorDim() const
  { 
    return this->ParamVectorDim();
  }
  
  /** Get parameter vector.
   */
  virtual void GetParamVector( CoordinateVector& v );

  /** Set parameter vector.
   */
  virtual void SetParamVector( CoordinateVector& v );

  /** Set parameter vector for a given transformation.
   */
  virtual void SetParamVector( CoordinateVector& v, const size_t xformIdx );

  /** Set single parameter value.
   */
  virtual void SetParameter( const size_t param, const Types::Coordinate value );

  /** Set single parameter value with separate xform and parameter index.
   */
  virtual void SetParameter( const size_t xform, const size_t param, const Types::Coordinate value );

  /** Evaluate functional with given parameter vector.
   * This function sets the current parameter vector, reformats all image data
   * into template space according to the current transformations, and calls
   * Evaluate() to compute the entropy-based functional value.
   *\param v Parameter vector.
   *\return Const function value for given parameters.
   */
  virtual Self::ReturnType EvaluateAt ( CoordinateVector& v );

  /** Compute functional value and gradient.
   *\param v Parameter vector.
   *\param g The extimated gradient of the functional is stored in this vector.
   *\param step Step size for finite difference gradient approximation. Default
   *  is 1 mm.
   *\return Const function value for given parameters.
   */
  virtual Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 );

  /** Allocate storage for reformatted images etc.
   * This function must be called AFTER setting template grid and target images, but BEFORE
   * any calls to Evaluate, EvaluateAt, or EvaluateWithGradient.
   */
  virtual void AllocateStorage();

  /// Write all images for debug purposes.
  void DebugWriteImages();

protected:
  /// Flag for freeing and re-reading original images if using smoothed data.
  bool m_FreeAndRereadImages;

  /// Flag for enforcing zero-sum parameter changes.
  bool m_ForceZeroSum;
  
  /// Restrict zero-sum computation to first N images.
  size_t m_ForceZeroSumFirstN;

  /// Currently active images from index.
  size_t m_ActiveImagesFrom;

  /// Currently active images to index (plus 1).
  size_t m_ActiveImagesTo;

  /// Currently active transformations from index.
  size_t m_ActiveXformsFrom;

  /// Currently active transformations to index (plus 1).
  size_t m_ActiveXformsTo;

  /// Enforce gradient to be zero-sum over all images.
  virtual void ForceZeroSumGradient( CoordinateVector& g ) const;

  /// Number of pixels in template.
  size_t m_TemplateNumberOfPixels;

  /// Number of samples drawn from the pixels in template.
  size_t m_TemplateNumberOfSamples;

  /// Template grid (not pixel data).
  UniformVolume::SmartPtr m_TemplateGrid;

  /// Flag for use of template pixel data in registration.
  bool m_UseTemplateData;

  /// Prepared (smoothed, scaled etc.) data of the template image if used in registration.
  std::vector<byte> m_TemplateData;

  /// Vector of image volumes with pre-scaled pixel values.
  std::vector<UniformVolume::SmartPtr> m_ImageVector;

  /// Vector of original image volumes.
  std::vector<UniformVolume::SmartPtr> m_OriginalImageVector;

  /// Vector of transformations
  std::vector<Xform::SmartPtr> m_XformVector;

  /// Probabilistic sample count.
  float m_ProbabilisticSampleDensity;

  /// Pixel indices of probabilistic samples.
  std::vector<size_t> m_ProbabilisticSamples;

  /** Number of iterations (calls to Evaluate()) after which probabilistic 
   * samples are updated.
   */
  int m_ProbabilisticSampleUpdatesAfter;

  /** Current number of iterations since last update of probabilistic samples.
   */
  int m_ProbabilisticSampleUpdatesSince;

  /** Update probablistic samples.
   * This function generates a new list of random pixel indices for sampling
   * the images. It is called every m_ProbabilisticSampleUpdatesAfter calls to
   * Evaluate().
   */
  virtual void UpdateProbabilisticSamples();

  /** Interpolate all moving images.
   * By default, this only calls InterpolateImage() for each image. However,
   * this functionality can be overriden for better efficiency, for example
   * in the distributed MPI implementation.
   */
  virtual void InterpolateAllImages();

  /** Interpolate given moving image to template.
   *\param idx Index of of to reformat to template. This also determines which
   *  transformation is used.
   *\param destination The reformatted pixel data is stored in this array.
   *  Sufficient memory (for as many pixels as there are in the template grid)
   *  must be allocated there.
   */
  virtual void InterpolateImage( const size_t idx, byte* const destination ) = 0;

  /// Vector of reformatted and rescaled image data.
  std::vector<byte*> m_Data;

  /// Temporary data allocated at correct size of template grid.
  std::vector<byte> m_TempData;

  /// Kernel width in [mm] for Gaussian smoothing of target images.
  Types::Coordinate m_GaussianSmoothImagesSigma;

  /// Value used to mark regions outside the FOV.
  static const byte m_PaddingValue = 255;

  /// User-defined value to fill regions outside FOV.
  byte m_UserBackgroundValue;
  
  /// Flag for user-defined background value.
  bool m_UserBackgroundFlag;

  /// Number of parameters per transformation..
  size_t m_ParametersPerXform;

  /// Repeat histogram-based intensity matching after each stage.
  bool m_RepeatIntensityHistogramMatching;

  /// Update probabilistic sample table..
  virtual bool Wiggle();

#ifdef CMTK_BUILD_MPI
  /// MPI process rank.
  int m_RankMPI;

  /// MPI process count.
  int m_SizeMPI;
#endif

  /// Prepare data for one image.
  virtual UniformVolume* PrepareSingleImage( UniformVolume::SmartPtr& image );

  /// Smooth and pre-scale target images.
  virtual void PrepareTargetImages();

private:
  /// Copy template data from TypedArray to byte vector.
  void CopyTemplateData();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkGroupwiseRegistrationFunctionalBase_h_included_
