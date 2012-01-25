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

#ifndef __cmtkImagePairNonrigidRegistrationFunctional_h_included_
#define __cmtkImagePairNonrigidRegistrationFunctional_h_included_

#include <cmtkconfig.h>

#include <Registration/cmtkImagePairRegistrationFunctional.h>
#include <Registration/cmtkImagePairSimilarityMeasure.h>

#include <System/cmtkThreads.h>
#include <System/cmtkThreadPool.h>

#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkJointHistogram.h>

#ifdef HAVE_IEEEFP_H
#  include <ieeefp.h>
#endif

#include <math.h>

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

  /// Virtual destructor.
  virtual ~ImagePairNonrigidRegistrationFunctional ();

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

  /** Active coordinate directions.
   */
  cmtkGetSetMacroString(ActiveCoordinates);

  /** Weight of the Jacobian constraint relative to voxel similarity measure.
   * If this is zero, only the voxel-based similarity will be computed.
   */
  cmtkGetSetMacroDefault(double,JacobianConstraintWeight,0);

  /** Weight of the grid energy relative to voxel similarity measure.
   * If this is zero, only the voxel-based similarity will be computed. If
   * equal to one, only the grid energy will be computed.
   */
  cmtkGetSetMacroDefault(double,GridEnergyWeight,0);

  /** Set Warp transformation.
   * This virtual function will be overridden by the derived classes that add
   * the actual warp transformation as a template parameters. It serves as a
   * common access point to update the warp transformation after construction
   * of the functional.
   */
  virtual void SetWarpXform( SplineWarpXform::SmartPtr& warp ) = 0;

  /// Set inverse transformation.
  void SetInverseTransformation( SplineWarpXform::SmartPtr& inverseTransformation ) 
  {
    this->m_InverseTransformation = inverseTransformation;
  }

  /// Set inverse consistency weight
  void SetInverseConsistencyWeight( const double inverseConsistencyWeight ) 
  {
    this->m_InverseConsistencyWeight = inverseConsistencyWeight;
  }

  /// Get parameter stepping in milimeters.
  virtual Types::Coordinate GetParamStep( const size_t idx, const Types::Coordinate mmStep = 1 ) const 
  {
    return this->m_Warp->GetParamStep( idx, this->m_FloatingSize, mmStep );
  }

  /// Return the transformation's parameter vector dimension.
  virtual size_t ParamVectorDim() const 
  {
    return this->m_Warp->ParamVectorDim();
  }

  /// Return the number of variable parameters of the transformation.
  virtual size_t VariableParamVectorDim() const 
  {
    return this->m_Warp->VariableParamVectorDim();
  }

  /// Return parameter vector.
  virtual void GetParamVector ( CoordinateVector& v ) 
  {
    this->m_Warp->GetParamVector( v );
  }

  /// Constructor function.
  static ImagePairNonrigidRegistrationFunctional* Create
  ( const int metric, UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume, const Interpolators::InterpolationEnum interpolation );
  
protected:
  /// Constructor.
  ImagePairNonrigidRegistrationFunctional( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating );

  /// Array of warp transformation objects for the parallel threads.
  std::vector<SplineWarpXform::SmartPtr> m_ThreadWarp;

  /// Array of storage for simultaneously retrieving multiple deformed vectors.
  Vector3D **m_ThreadVectorCache;

  /** Number of actual parallel threads used for computations.
   * All duplicated data structures are generated with the multiplicity given
   * by this value. It is determined from Threads when the object is first
   * instanced. It cannot be changed afterwards.
   */
  size_t m_NumberOfThreads;

  /// Number of parallel tasks.
  size_t m_NumberOfTasks;

  /// Baseline transformed volume.
  Types::DataItem *m_WarpedVolume;

  /// Shortcut variables for x, y, z dimension of the reference image.
  DataGrid::IndexType::ValueType m_DimsX, m_DimsY, m_DimsZ;
  
  /// Shorcut variables for x and y dimensions of the floating image.
  DataGrid::IndexType::ValueType m_FltDimsX, m_FltDimsY;

  /// Pointer to the local warp transformation.
  SplineWarpXform::SmartPtr m_Warp;

  /// Optional inverse transformation for inverse-consistent deformation.
  SplineWarpXform::SmartPtr m_InverseTransformation;

  /// Weight for inverse consistency constraint.
  double m_InverseConsistencyWeight;

  /// Return weighted combination of voxel similarity and grid energy.
  Self::ReturnType WeightedTotal( const Self::ReturnType metric, const SplineWarpXform& warp ) const 
  {
    double result = metric;
    if ( this->m_JacobianConstraintWeight > 0 ) 
      {
      result -= this->m_JacobianConstraintWeight * warp.GetJacobianConstraint();
      } 
    
    if ( this->m_GridEnergyWeight > 0 ) 
      {
      result -= this->m_GridEnergyWeight * warp.GetGridEnergy();
      }
    
    if ( !finite( result ) ) 
      return -FLT_MAX;
    
    if ( this->m_MatchedLandmarkList ) 
      {
      result -= this->m_LandmarkErrorWeight * warp.GetLandmarksMSD( this->m_MatchedLandmarkList );
      }

    if ( this->m_InverseTransformation ) 
      {
      result -= this->m_InverseConsistencyWeight * warp.GetInverseConsistencyError( this->m_InverseTransformation, this->m_ReferenceGrid );
      }
    
    return static_cast<Self::ReturnType>( result );
  }
  
  /// Return weighted combination of similarity and grid energy derivatives.
  void WeightedDerivative( double& lower, double& upper, SplineWarpXform& warp, const int param, const Types::Coordinate step ) const;

  /** Regularize the deformation.
   */
  cmtkGetSetMacroDefault(bool,Regularize,false);

  /// Dimension of warp parameter vector
  size_t Dim;

  /** Parameter scaling vector.
   * This array holds the scaling factors for all warp parameters as returned
   * by the transformation class. These factors can be used to equalized all
   * parameter modifications during gradient computation etc.
   */
  std::vector<Types::Coordinate> m_StepScaleVector;

  /** Volume of influence table.
   * This array holds the precomputed volumes of influence for all 
   * transformation parameters. Six successive numbers per parameter define the
   * voxel range with respect to the reference colume grid that is affected by
   * the respective parameter.
   */
  DataGrid::RegionType *VolumeOfInfluence;

  /// Coordinate domain of the reference image.
  UniformVolume::CoordinateRegionType m_ReferenceDomain;

  /// Make smart pointer class friend so we can keep destructor protected.
  friend class SmartPointer<Self>;
};

//@}

} // namespace cmtk

#endif // __cmtkImagePairNonrigidRegistrationFunctional_h_included_
