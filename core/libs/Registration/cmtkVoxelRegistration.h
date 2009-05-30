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

#ifndef __cmtkVoxelRegistration_h_included_
#define __cmtkVoxelRegistration_h_included_

#include <cmtkconfig.h>

#include <cmtkMacros.h>
#include <cmtkTypes.h>

#include <cmtkXform.h>
#include <cmtkAffineXform.h>
#include <cmtkFunctional.h>
#include <cmtkUniformVolume.h>
#include <cmtkMatchedLandmarkList.h>
#include <cmtkVector.h>

#include <cmtkRegistrationCallback.h>
#include <cmtkOptimizer.h>

#include <cmtkTimers.h>

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
  igsGetSetMacro(int,Metric);

  /// Optimization algorithm to use.
  igsGetSetMacro(int,Algorithm);

  /// Exploration, i.e. initial step size.
  igsGetSetMacro(double,Exploration);

  /// Accuracy, i.e. final step size.
  igsGetSetMacro(double,Accuracy);

  /** Coarsest resolution to resample image data to.
   * If this value is unset, ie. less than or equal to zero, then the coarsest
   * resolution is automatically computed from the initial step size
   * (Exploration).
   */
  double CoarsestResolution;

  /// Flag whether the last resolution level uses the original images.
  igsGetSetMacro(bool,UseOriginalData);
 
 /// Factor between optimization step sizes.
  double OptimizerStepFactor;

  /// Use maximum norm instead of Euclid where applicable.
  bool UseMaxNorm;

  /// Sampling, i.e. last non-original resolution.
  igsGetSetMacro(Types::Coordinate,Sampling);

  /// Name of protocol file.
  igsGetSetMacroString(Protocol);

  /// First data volume.
  igsGetSetMacro(UniformVolume::SmartPtr,Volume_1);

  /// Data class of first data volume.
  igsGetSetMacro(DataClass,DataClass_1);

  /// Second data volume.
  igsGetSetMacro(UniformVolume::SmartPtr,Volume_2);

  /// Data class of second data volume.
  igsGetSetMacro(DataClass,DataClass_2);

  /// Lower threshold flag for image 1.
  bool m_ThreshMin1;

  /// Lower threshold value for image 1.
  float m_ThreshMinValue1;

  /// Lower threshold flag for image 2.
  bool m_ThreshMin2;

  /// Lower threshold value for image 2.
  float m_ThreshMinValue2;

  /// Upper threshold flag for image 1.
  bool m_ThreshMax1;

  /// Upper threshold value for image 1.
  float m_ThreshMaxValue1;

  /// Upper threshold flag for image 2.
  bool m_ThreshMax2;

  /// Upper threshold value for image 2.
  float m_ThreshMaxValue2;

  /// Weighting factor of landmark registration error vs. image similarity.
  igsGetSetMacro(float,LandmarkErrorWeight);

  /// Matched landmarks list.
  igsGetSetMacro(MatchedLandmarkList::SmartPtr,LandmarkList);

  /// Flag whether model and reference are exchanged.
  bool SwitchVolumes;

  /// Pointer to callback object.
  igsGetSetMacro(RegistrationCallback::SmartPtr,Callback);

  /// Initial transformation.
  igsGetSetMacro(AffineXform::SmartPtr,InitialXform);

  /// FLag whether initial transformation is an inverse.
  igsGetSetMacro(bool,InitialXformIsInverse);

  /// Current / final transformation.
  Xform::SmartPtr Xform;

  /// Stack of functional objects for the resolution steps.
  std::stack<Functional::SmartPtr> FunctionalStack;

  /// Pointer to optimizer object.
  Optimizer::SmartPtr Optimizer;

  /**@name Member functions to be overwritten.
   */
  //@{
  /** Initialize registration.
   * This function is called by Register before any other operations. It can
   * be overloaded to open status dialog windows, etc. Derived implementations
   * should call their base class' InitRegistration first.
   *@return 
   */
  virtual CallbackResult InitRegistration ();

  /** Report registration progress.
   * This function is called repeatedly during creation of the resampled
   * data volumes. It is also called before optimizing every resolution
   * level. It is not, however, called DURING optimization. For that, use
   * a suitable callback object.
   *@return If this function returns 1, registration is terminated as soon
   * as possible.
   *@see CreateCallback
   */
  virtual CallbackResult ReportProgress ( const char*, const int ) const 
  { return CALLBACK_OK; };

  /** Output registration result.
   * This function is called after finishing registration. It can overloaded
   * to report the resulting transformation, export it to an encapsulating
   * application, etc...
   *@param v The vector of resulting transformation parameters.
   */
  virtual void OutputResult ( const CoordinateVector* ) const {};
  
  /** Finalize registration.
   * This function is called after registration has been terminated. It can
   * be used to destroy progress dialog windows, free memory etc. Its last
   * operation should be a call to the respective parent class' implementation.
   */
  virtual void DoneRegistration ( const CoordinateVector* v = NULL);

  /** Enter resolution level.
   * This function is called before entering each resolution level. It can
   * be used to update status displays etc.
   *@param index Index of the current resolution level. 0 is first (coarsest),
   * subsequent (finer) resolutions have increasing numbers.
   *@param total Total number of resolution levels.
   */
  virtual void EnterResolution( CoordinateVector::SmartPtr&, Functional::SmartPtr&, const int, const int );

  /** Finish resolution level.
   * This function is called after every resolution level. It should do any
   * necessary cleanups resulting from the previous call to EnterRegistration.
   *@return If the current level is finished, 1 is returned. Otherwise, ie.
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
   *@return 1 if registration was terminated by a user interrupt, 0 if
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
