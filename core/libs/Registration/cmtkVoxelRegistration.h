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

#ifndef __cmtkVoxelRegistration_h_included_
#define __cmtkVoxelRegistration_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkMacros.h>
#include <Base/cmtkTypes.h>

#include <Base/cmtkXform.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkFunctional.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkVector.h>

#include <Registration/cmtkRegistrationCallback.h>
#include <Registration/cmtkOptimizer.h>
#include <Base/cmtkInterpolator.h>

#include <System/cmtkTimers.h>
#include <System/cmtkCommandLine.h>

#include <IO/cmtkClassStreamOutput.h>

#include <stack>
#include <string.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Generic multiresolution voxel-registration class.
 * By implementing member functions to retrieve parameters and report results
 * in derived classes, registration can be integrated into various 
 * environments.
 */
class VoxelRegistration 
{
protected:
  /// Metric to use.
  cmtkGetSetMacro(int,Metric);

  /// Optimization algorithm to use.
  cmtkGetSetMacro(int,Algorithm);

  /// Exploration, i.e. initial step size.
  cmtkGetSetMacro(double,Exploration);

  /// Accuracy, i.e. final step size.
  cmtkGetSetMacro(double,Accuracy);

  /** Coarsest resolution to resample image data to.
   * If this value is unset, ie. less than or equal to zero, then the coarsest
   * resolution is automatically computed from the initial step size
   * (Exploration).
   */
  double CoarsestResolution;

  /// Flag whether the last resolution level uses the original images.
  cmtkGetSetMacro(bool,UseOriginalData);
 
 /// Factor between optimization step sizes.
  double OptimizerStepFactor;

  /// Use maximum norm instead of Euclid where applicable.
  bool UseMaxNorm;

  /// Threshold for terminating optimization based on changes of the target function.
  Optimizer::ReturnType m_DeltaFThreshold;

  /// Sampling, i.e. last non-original resolution.
  cmtkGetSetMacro(Types::Coordinate,Sampling);

  /// Name of protocol file.
  cmtkGetSetMacroString(Protocol);

  /// First data volume.
  cmtkGetSetMacro(UniformVolume::SmartPtr,Volume_1);

  /// Second data volume.
  cmtkGetSetMacro(UniformVolume::SmartPtr,Volume_2);

  /** Reference data volume.
   * This is a pointer to the actual reference volume, which is either Volume_1 or Volume_2 above,
   * depending on whether registration was instructed to switch the two or not.
   */
  cmtkGetSetMacro(UniformVolume::SmartPtr,ReferenceVolume);

  /** Floating data volume.
   * This is a pointer to the actual floating volume, which is either Volume_2 or Volume_1 above,
   * depending on whether registration was instructed to switch the two or not.
   */
  cmtkGetSetMacro(UniformVolume::SmartPtr,FloatingVolume);

  /// Local class for preprocessing image data, e.g., by histogram operations, thresholding, and cropping.
  class ImagePreprocessor
  {
  public:
    /// Data class string ("grey", "labels", or "binary")
    const char* m_DataClassString;

    /// Data class ID.
    DataClass m_DataClass;
    
    /// Flag for pixel padding.
    bool m_PaddingFlag;
    
    /// Padding value.
    Types::DataItem m_PaddingValue;

    /// Lower threshold flag.
    bool m_LowerThresholdActive;

    /// Lower threshold value.
    Types::DataItem m_LowerThresholdValue;
  
    /// Upper threshold flag.
    bool m_UpperThresholdActive;

    /// Upper threshold value.
    Types::DataItem m_UpperThresholdValue;

    /// Flag for image histogram pruning.
    bool m_UsePruneHistogramBins;

    /// Prune histogram for image: number of target bins.
    unsigned int m_PruneHistogramBins;

    /// Flag for histogram equalization.
    bool m_HistogramEqualization;

    /// Flag for application of Sobel edge detection filter.
    bool m_SobelFilter;

    /// Crop region in index coordinates.
    const char* m_CropIndex;

    /// Crop region in world coordinates.
    const char* m_CropWorld;

    /// Flag for auto cropping.
    bool m_AutoCropFlag;

    /// Auto cropping level.
    Types::DataItem m_AutoCropLevel;
    
    /// Constructor.
    ImagePreprocessor( const char* name /*!< There are two preprocessors, for reference and floating image: this parameter names a parameter group for this instance.*/,
		       const char* key /*!< This parameter gives a string key that is appended to each command line option so that reference and floating preprocessors do not collide.*/ );
    
    /// Attach this preprocessor to a command line parse.
    void AttachToCommandLine( CommandLine& cl /*!< The command line object to add our options to.*/ ); 
    
    /// Get pre-processed image from original image.
    UniformVolume::SmartPtr GetProcessedImage( const UniformVolume* original );
    
    /// Write settings of this object to class stream for archiving.
    void WriteSettings( ClassStreamOutput& stream ) const;

  private:
    /// Store the name that identifies this instance ("Reference" or "Floating")
    const char* m_Name;

    /// Store the key that identifies this instance ("ref" or "flt")
    const char* m_Key;
  };

  /// Image preprocessor for reference image.
  ImagePreprocessor m_PreprocessorRef;

  /// Image preprocessor for floating image.
  ImagePreprocessor m_PreprocessorFlt;
  
  /// Flag whether model and reference are exchanged.
  bool SwitchVolumes;

  /// Pointer to callback object.
  cmtkGetSetMacro(RegistrationCallback::SmartPtr,Callback);

  /// Initial transformation.
  cmtkGetSetMacro(AffineXform::SmartPtr,InitialTransformation);

  /// FLag whether initial transformation is an inverse.
  cmtkGetSetMacro(bool,InitialXformIsInverse);

  /// Current / final transformation.
  Xform::SmartPtr m_Xform;

  /// Stack of functional objects for the resolution steps.
  std::stack<Functional::SmartPtr> FunctionalStack;

  /// Pointer to optimizer object.
  Optimizer::SmartPtr m_Optimizer;

  /**\name Member functions to be overwritten.
   */
  //@{
  /** Initialize registration.
   * This function is called by Register before any other operations. It can
   * be overloaded to open status dialog windows, etc. Derived implementations
   * should call their base class' InitRegistration first.
   *\return 
   */
  virtual CallbackResult InitRegistration ();

  /** Output registration result.
   * This function is called after finishing registration. It can overloaded
   * to report the resulting transformation, export it to an encapsulating
   * application, etc...
   */
  virtual void OutputResult ( const CoordinateVector* /*!< The vector of resulting transformation parameters. */, 
			      const CallbackResult = CALLBACK_OK /*!< The interrupt status - this allows the output function to determine whether computation finished or was interrupted. */ ) {}
  
  /** Finalize registration.
   * This function is called after registration has been terminated. It can
   * be used to destroy progress dialog windows, free memory etc. Its last
   * operation should be a call to the respective parent class' implementation.
   */
  virtual void DoneRegistration ( const CoordinateVector* v = NULL);

  /** Enter resolution level.
   * This function is called before entering each resolution level. It can
   * be used to update status displays etc.
   *\param v Current parameter vector.
   *\param f Functional for next level.
   *\param idx Index of the current resolution level. 0 is first (coarsest),
   * subsequent (finer) resolutions have increasing numbers.
   *\param total Total number of resolution levels.
   */
  virtual void EnterResolution( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, const int idx, const int total );

  /** Finish resolution level.
   * This function is called after every resolution level. It should do any
   * necessary cleanups resulting from the previous call to EnterRegistration.
   *\return If the current level is finished, 1 is returned. Otherwise, ie.
   * if the derived class requests another run of the same level, 0 may be
   * returned. This is used for example by the affine registration in order
   * to make repeated runs of the same level with different numbers of degrees
   * of freedom. Be careful not to create any inifinite loops.
   */
  virtual int DoneResolution( CoordinateVector::SmartPtr&, Functional::SmartPtr&, const int, const int ) { return 1; }
  //@}

public:
  /// Exception class.
  class ConstructorFailed {};

  /** Default constructor. 
   */
  VoxelRegistration ();

  /** Destructor.
   */
  virtual ~VoxelRegistration ();

  /** Do registration.
   * This function must be called to start the multiresolution optimization
   * using the specified parameters.
   *\return 1 if registration was terminated by a user interrupt, 0 if
   * registration was finished.
   */
  virtual CallbackResult Register ();

  /** Return total elapsed process time.
   */
  double GetTotalElapsedTime() const 
  {
    return cmtk::Timers::GetTimeProcess() - TimeStartRegistration;
  }
  
  /** Return elapsed process time during current level.
   */
  double GetLevelElapsedTime() const 
  {
    return cmtk::Timers::GetTimeProcess() - TimeStartLevel;
  }
  
  /** Return total elapsed walltime.
   */
  double GetTotalElapsedWalltime() const 
  {
    return cmtk::Timers::GetWalltime() - WalltimeStartRegistration;
  }
  
  /** Return elapsed walltime during current level.
   */
  double GetLevelElapsedWalltime() const 
  {
    return cmtk::Timers::GetWalltime() - WalltimeStartLevel;
  }
  
  /** Return total elapsed thread time.
   */
  double GetThreadTotalElapsedTime() const 
  {
    return cmtk::Timers::GetTimeThread() - ThreadTimeStartRegistration;
  }
  
  /** Return elapsed thread time during current level.
   */
  double GetThreadLevelElapsedTime() const 
  {
    return cmtk::Timers::GetTimeThread() - ThreadTimeStartLevel;
  }

private:
  /** Time of registration start.
   * This is used as the reference for absolute computation time calculation.
   */
  double TimeStartRegistration;

  /** Time of entering the current resolution level.
   * This is used as the reference for per-level computation time calculation.
   */
  double TimeStartLevel;

  /** Reference walltime of registration start.
   * This is used as the reference for absolute computation time calculation.
   */
  double WalltimeStartRegistration;

  /** Reference walltime of entering the current resolution level.
   * This is used as the reference for per-level computation time calculation.
   */
  double WalltimeStartLevel;

  /** Time of registration start.
   * This is used as the reference for absolute computation time calculation.
   */
  double ThreadTimeStartRegistration;

  /** Time of entering the current resolution level.
   * This is used as the reference for per-level computation time calculation.
   */
  double ThreadTimeStartLevel;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVoxelRegistration_h_included_
