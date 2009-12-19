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

#ifndef __cmtkImagePairNonrigidRegistrationFunctional_h_included_
#define __cmtkImagePairNonrigidRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <cmtkImagePairRegistrationFunctional.h>

#include <cmtkThreads.h>
#include <cmtkThreadPool.h>

#include <cmtkWarpXform.h>
#include <cmtkJointHistogram.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Common base class for all elastic registration functionals.
 * This class holds all members that are not related to the effective metric
 * and therefore need not be present in the derived template class.
 */
class ImagePairNonrigidRegistrationFunctional : 
  /// Inherit basic image pair registration functional class.
  public ImagePairRegistrationFunctional 
{
public:
  /// This class.
  typedef ImagePairNonrigidRegistrationFunctional Self;

  /// Smart pointer to this class.
  typedef SmartPointer<Self> SmartPtr;

  /// Superclass.
  typedef ImagePairRegistrationFunctional Superclass;

  /** Set Warp transformation.
   * This virtual function will be overridden by the derived classes that add
   * the actual warp transformation as a template parameters. It serves as a
   * common access point to update the warp transformation after construction
   * of the functional.
   */
  virtual void SetWarpXform( WarpXform::SmartPtr& warp ) = 0;

  /// Set flag and value for forcing pixels outside the floating image.
  virtual void SetForceOutside( const bool flag = true, const Types::DataItem value = 0 ) = 0;

  /// Destructor.
  virtual ~ImagePairNonrigidRegistrationFunctional ();

  /// Constructor function.
  static ImagePairNonrigidRegistrationFunctional* Create
  ( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume, const Interpolators::InterpolationEnum interpolation );
  
protected:
  /// Constructor.
  ImagePairNonrigidRegistrationFunctional( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating );

  /** Set active and passive warp parameters adaptively.
   * If this flag is set, the functional will adaptively determine active and
   * passive parameters of the warp transformation prior to gradient 
   * computation.
   */
  cmtkGetSetMacroDefault(bool,AdaptiveFixParameters,true);

  /** Set threshold factor for selecting passive warp parameters adaptively.
   * If the flag AdaptiveFixParameters is set, this value determines the
   * threshold by which active vs. passive parameters are selected. All
   * control points are set to passive for which the local region entropy is
   * below this factor times sum of min and max region entropy. The default
   * value is 0.5.
   */
  cmtkGetSetMacro(double,AdaptiveFixThreshFactor);

  /** Weight of the Jacobian constraint relative to voxel similarity measure.
   * If this is zero, only the voxel-based similarity will be computed.
   */
  cmtkGetSetMacroDefault(double,JacobianConstraintWeight,0);

  /** Spatial map of relative (tissue-specific) incompressibility constraint.
   */
  cmtkGetSetMacro(DataGrid::SmartPtr,IncompressibilityMap);

  /** Weight of the grid energy relative to voxel similarity measure.
   * If this is zero, only the voxel-based similarity will be computed. If
   * equal to one, only the grid energy will be computed.
   */
  cmtkGetSetMacroDefault(double,GridEnergyWeight,0);

  /** Regularize the deformation.
   */
  cmtkGetSetMacroDefault(bool,Regularize,false);

  /** Warp's fixed parameters need to be updated.
   * This flag is set when the warp transformation is set or modified. It
   * signals that the active and passive parameters of the transformation
   * will have to be updated before the next gradient computation.
   */
  bool WarpNeedsFixUpdate;

  /// Histogram used for consistency computation.
  JointHistogram<unsigned int>::SmartPtr m_ConsistencyHistogram;

  /// Dimension of warp parameter vector
  size_t Dim;

  /** Parameter scaling vector.
   * This array holds the scaling factors for all warp parameters as returned
   * by the transformation class. These factors can be used to equalized all
   * parameter modifications during gradient computation etc.
   */
  Types::Coordinate *StepScaleVector;

  /** Volume of influence table.
   * This array holds the precomputed volumes of influence for all 
   * transformation parameters. Six successive numbers per parameter define the
   * voxel range with respect to the reference colume grid that is affected by
   * the respective parameter.
   */
  Rect3D *VolumeOfInfluence;

  /// Coordinate of the beginning of the reference colume crop area.
  Vector3D ReferenceFrom;

  /// Coordinate of the end of the reference colume crop area.
  Vector3D ReferenceTo;

  /// Storage for simultaneously retrieving multiple deformed vectors.
  Vector3D *VectorCache;
};

//@}

} // namespace cmtk

#endif // __cmtkImagePairNonrigidRegistrationFunctional_h_included_
