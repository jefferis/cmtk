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

#ifndef __cmtkSplineWarpGroupwiseRegistrationRMIFunctional_h_included_
#define __cmtkSplineWarpGroupwiseRegistrationRMIFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkGroupwiseRegistrationRMIFunctional.h>

#include <cmtkSmartPtr.h>
#include <cmtkThreads.h>

#include <cmtkUniformVolume.h>
#include <cmtkSplineWarpXform.h>
#include <cmtkHistogram.h>

#include <vector>

#include <cmtkClassStream.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Functional for spline warp groupwise registration.
 */
class SplineWarpGroupwiseRegistrationRMIFunctional : 
  public GroupwiseRegistrationRMIFunctional<SplineWarpXform>
{
public:
  /// Type of parent class.
  typedef GroupwiseRegistrationRMIFunctional<SplineWarpXform> Superclass;

  /// Type of this class.
  typedef SplineWarpGroupwiseRegistrationRMIFunctional Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  SplineWarpGroupwiseRegistrationRMIFunctional();

  /// Destructor.
  virtual ~SplineWarpGroupwiseRegistrationRMIFunctional();

  /** Initialize spline warp transformations.
   */
  void InitializeXforms( const Types::Coordinate gridSpacing, std::vector<AffineXform::SmartPtr> initialAffineXformsVector );

  /** Initialize spline warp transformations.
   */
  void InitializeXforms( const Types::Coordinate gridSpacing )
  {
    this->InitializeXforms( gridSpacing, this->m_InitialAffineXformsVector );
  }

  /// Refine transformation control point grids.
  void RefineTransformationGrids();

  /// Set partial gradient mode.
  void SetPartialGradientMode( const bool partialGradientMode = false, const Types::Coordinate partialGradientThreshold = 0.0 )
  {
    this->m_PartialGradientMode = partialGradientMode;
    this->m_PartialGradientThreshold = partialGradientThreshold;
  }

  /// Set deactivate uninformative control points mode.
  void SetDeactivateUninformativeMode( const bool dum = true )
  {
    this->m_DeactivateUninformativeMode = dum;
  }

  /** Set range of currently active transformations.
   * Call inherited function, then update local step size array.
   */
  virtual void SetActiveXformsFromTo( const size_t from, const size_t to )
  {
    this->Superclass::SetActiveXformsFromTo( from, to );
    this->UpdateParamStepArray();
  }

  /// Call inherited function and allocate local storage.
  virtual void SetTemplateGrid( UniformVolume::SmartPtr& templateGrid, const int downsample = 1,const bool useTemplateData = false );
    
  /// Evaluate functional with currently set parameters.
  virtual Self::ReturnType Evaluate();

  /// Evaluate functional and set parameters.
  virtual Self::ReturnType EvaluateAt( CoordinateVector& v )
  {
    return this->Superclass::EvaluateAt( v );
  }

  /** Compute functional value and gradient.
   *\param v Parameter vector.
   *\param g The extimated gradient of the functional is stored in this vector.
   *\param step Step size for finite difference gradient approximation. Default
   *  is 1 mm.
   *\return Const function value for given parameters.
   */
  virtual Self::ReturnType EvaluateWithGradient( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step = 1 );

  /// Reset gradient dimensions, etc.
  virtual bool Wiggle();

protected:
  /** Interpolate given moving image to template.
   * This function overrides the interpolation function provided by the base
   * class. It makes use of the fact that affine transformations preserve
   * parallel lines for more efficient computation.
   *\param idx Index of of to reformat to template. This also determines which
   *  transformation is used.
   *\param destination The reformatted pixel data is stored in this array.
   *  Sufficient memory (for as many pixels as there are in the template grid)
   *  must be allocated there.
   */
  virtual void InterpolateImage( const size_t idx, byte* const destination );

private:
  /// Flag for fast warp mode, i.e., reduced control point influence volume.
  bool m_WarpFastMode;

  /** Flag for partial gradient computation.
   * If this is set, gradient components under a given threshold are deactivated
   * and not used for gradient approximation.
   */
  bool m_PartialGradientMode;

  /// Threshold for partial gradient computation.
  Types::Coordinate m_PartialGradientThreshold;

  /// Deactivate uninformative control points mode.
  bool m_DeactivateUninformativeMode;

  /// List of flags for deactivated control points.
  std::vector<bool> m_ActiveControlPointFlags;

  /// Update volumes of influence for warp parameters.
  virtual void UpdateVolumesOfInfluence();

  /// Number of deactivated control points.
  size_t m_NumberOfActiveControlPoints;

  /// Update deactivated control points.
  virtual void UpdateActiveControlPoints();

  /// Local information measure for neighborhood of each control point.
  std::vector<byte> m_InformationByControlPoint;

  /// Flag whether information by control point needs to be updated.
  bool m_NeedsUpdateInformationByControlPoint;

  /// Update local information by control point.
  virtual void UpdateInformationByControlPoint();

  /// Update control point schedule for gradient approximation.
  virtual void UpdateControlPointSchedule();

  /** Update list of probabilistic samples.
   * Call inherited function, then determine which parameters of the current
   * warp affect samples in the list and deactivate all others.
   */
  virtual void UpdateProbabilisticSamples();

  /// Initial affine transformations.
  std::vector<AffineXform::SmartPtr> m_InitialAffineXformsVector;

  /// Current parameter steppings for the warp parameters.
  std::vector<Types::Coordinate> m_ParamStepArray;

  /// Update parameter steppings for the warp parameters.
  virtual bool UpdateParamStepArray();

  /// Volumes of influence for the warp parameters.
  std::vector<Rect3D> m_VolumeOfInfluenceArray;

  /// Thread function parameters for image interpolation.
  class InterpolateImageThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParameters<Self>
  {
  public:
    /// Index of the image to be interpolated.
    size_t m_Idx;

    /// Pointer to storage that will hold the reformatted pixel data.
    byte* m_Destination;
  };

  /// Image interpolation thread function.
  static CMTK_THREAD_RETURN_TYPE InterpolateImageThread( void* args );

  /// Processing schedule for overlap-free parallel processing of control points.
  std::vector<int> m_ControlPointSchedule;

  /// Maximum number of concurrent jobs working on warps that is still overlap-free.
  size_t m_ControlPointScheduleOverlapFreeMaxLength;

  /// Thread function parameters for image interpolation.
  class EvaluateLocalGradientThreadParameters : 
    /// Inherit from generic thread parameters.
    public ThreadParameters<Self>
  {
  public:
    /// Unique thread storage index.
    size_t m_ThreadStorageIndex;

    /// Global parameter step scale factor.
    Types::Coordinate m_Step;

    /// Pointer to array with computed local gradient components.
    Types::Coordinate* m_Gradient;
    
    /// Current metric value.
    Self::ReturnType m_MetricBaseValue;

#ifdef CMTK_BUILD_MPI
    /// Index of first control point to be computed by all threads in this iteration.
    size_t m_FirstIndexToCompute;
#endif
  };

  /** Thread function: Compute local gradient of the cost function for gradient approximation.
   * This function takes into consideration that in a spline warp, each control point
   * effects only a local neighborhood. It also groups the parameters by control
   * point and works over all images and x,y,z to speed things up substantially.
   */
  static CMTK_THREAD_RETURN_TYPE EvaluateLocalGradientThreadFunc( void* args );

#ifdef CMTK_BUILD_MPI
  void ReorderGradientComponents( Types::Coordinate *const dst, const Types::Coordinate* src, const size_t fromCpIdx, const size_t toCpIdx );
#endif

  /// Function friends.
  friend ClassStream& operator<<( ClassStream& stream, const SplineWarpGroupwiseRegistrationRMIFunctional& func );
  friend ClassStream& operator>>( ClassStream& stream, SplineWarpGroupwiseRegistrationRMIFunctional& func );
};

/// Class stream write function.
ClassStream& operator<<( ClassStream& stream, const SplineWarpGroupwiseRegistrationRMIFunctional& func );

/// Class stream read function.
ClassStream& operator>>( ClassStream& stream, SplineWarpGroupwiseRegistrationRMIFunctional& func );

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSplineWarpGroupwiseRegistrationRMIFunctional_h_included_
