/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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
//  $Revision: 1768 $
//
//  $LastChangedDate: 2010-05-27 15:16:01 -0700 (Thu, 27 May 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkGroupwiseRegistrationFunctionalXformTemplate_SplineWarp_h_included_
#define __cmtkGroupwiseRegistrationFunctionalXformTemplate_SplineWarp_h_included_

#include <Base/cmtkSplineWarpXform.h>

namespace
cmtk
{

/** Template specialization for groupwise nonrigid registration functionals.
 * This class is the specialization of the generic transformation-dependent
 * functional class template, specialized for nonrigid (B-spline FFD) transformations.
 *
 * As such, this class provides functionality such as: initialization of
 * FFDs from affine transformations, grid refinement, and deformation constraints.
 */
template<>
class GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform> : 
  /** Inherit from generic groupwise functional. */
  public GroupwiseRegistrationFunctionalXformTemplateBase<SplineWarpXform>
{
public:
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalXformTemplateBase<SplineWarpXform> Superclass;
  
  /// Type of this class.
  typedef GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform> Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /// Constructor.
  GroupwiseRegistrationFunctionalXformTemplate();

  /// Destructor.
  virtual ~GroupwiseRegistrationFunctionalXformTemplate() {};

  /// Initialize nonrigid from affine transformations.
  virtual void InitializeXformsFromAffine( const Types::Coordinate gridSpacing /*!< Control point grid spacing in real-world units*/,
					   std::vector<AffineXform::SmartPtr> initialAffineXformsVector /*!< Vector of initial affine coordinate transformations*/,
					   const bool exactSpacing = true /*!< If set, the control point spacing will be exactly as given in the first parameter*/ );
  
  /** Initialize spline warp transformations.
   */
  virtual void InitializeXforms( const Types::Coordinate gridSpacing /*!< Control point grid spacing in real-world units*/,
				 const bool exactSpacing = true /*!< If set, the control point spacing will be exactly as given in the first parameter*/ )
  {
    this->InitializeXformsFromAffine( gridSpacing, this->m_InitialAffineXformsVector, exactSpacing );
  }
  
  /// Refine transformation control point grids.
  virtual void RefineTransformationGrids();

  /// Set flag for exclusion of affine components in unbiased groupwise deformation.
  void SetForceZeroSumNoAffine( const bool noaffine = true )
  {
    this->m_ForceZeroSumNoAffine = noaffine;
  }

  /// Set partial gradient mode.
  void SetPartialGradientMode( const bool partialGradientMode = false, const float partialGradientThreshold = 0.0 )
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
  virtual void SetTemplateGrid( UniformVolume::SmartPtr& templateGrid, const int downsample = 1, const bool useTemplateData = false );
    
protected:
  /// Maximum number of pixels in any VOI.
  size_t m_MaximumNumberOfPixelsVOI;
  
  /// Maximum number of pixels per line in any VOI.
  size_t m_MaximumNumberOfPixelsPerLineVOI;

  /// Update volumes of influence for warp parameters.
  virtual void UpdateVolumesOfInfluence();

  /// Update deactivated control points.
  virtual void UpdateActiveControlPoints() = 0;
  
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

  /** Enforce gradient to be zero-sum over all images.
   * This function essentially calls the inherited function of the same name. However,
   * if this->m_ForceZeroSumNoAffine is true, then the initial affine transformations
   * of each warp are eliminated from the gradient prior to calling the inherited
   * function, and they are re-applied afterwards. This way, the unbiased property of
   * the transformation set is made invariant under the affine transformation components.
   */
  virtual void ForceZeroSumGradient( CoordinateVector& g ) const;

  /** "Wiggle" functional a little, i.e., by updating probabilistic sampling.
   *\return True if the functional actually changed slightly, false if no "wiggle"
   * has actually taken place. This signals the optimizer whether it's worth trying
   * the current stage again (true) or not (false).
   */
  bool Wiggle();

protected:
  /// Flag for correction of affine components in unbiased warp.
  bool m_ForceZeroSumNoAffine;

  /// Flag for fast warp mode, i.e., reduced control point influence volume.
  bool m_WarpFastMode;

  /// Weight for jacobian constraint term.
  float m_JacobianConstraintWeight;

  /// Weight for grid bending energy term.
  float m_BendingEnergyWeight;

  /** Flag for partial gradient computation.
   * If this is set, gradient components under a given threshold are deactivated
   * and not used for gradient approximation.
   */
  bool m_PartialGradientMode;

  /// Initial affine transformations.
  std::vector<AffineXform::SmartPtr> m_InitialAffineXformsVector;

  /// Rotation components of initial affine transformations.
  std::vector<AffineXform::SmartPtr> m_InitialRotationsVector;

  /// Current parameter steppings for the warp parameters.
  std::vector<Types::Coordinate> m_ParamStepArray;

  /// Update parameter steppings for the warp parameters.
  virtual bool UpdateParamStepArray();

  /// Volumes of influence for the warp parameters.
  std::vector<DataGrid::RegionType> m_VolumeOfInfluenceArray;

  /// Threshold for partial gradient computation.
  Types::Coordinate m_PartialGradientThreshold;

  /// Deactivate uninformative control points mode.
  bool m_DeactivateUninformativeMode;

  /// List of flags for deactivated control points.
  std::vector<bool> m_ActiveControlPointFlags;

  /// Number of deactivated control points.
  size_t m_NumberOfActiveControlPoints;

private:
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
  static void InterpolateImageThread( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t );

  friend ClassStream& operator<<( ClassStream& stream, const GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>& func );
  friend ClassStream& operator>>( ClassStream& stream, GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>& func );
};

ClassStream& operator<<( ClassStream& stream, const GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>& func );
ClassStream& operator>>( ClassStream& stream, GroupwiseRegistrationFunctionalXformTemplate<SplineWarpXform>& func );

} // namespace cmtk

#endif // #ifndef __cmtkGroupwiseRegistrationFunctionalXformTemplate_SplineWarp_h_included_
