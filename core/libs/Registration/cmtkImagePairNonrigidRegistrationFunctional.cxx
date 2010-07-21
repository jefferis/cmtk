/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include "cmtkImagePairNonrigidRegistrationFunctional.h"

#include "Base/cmtkSplineWarpXform.h"
#include "Base/cmtkInterpolator.h"

#include "Registration/cmtkImagePairNonrigidRegistrationFunctionalTemplate.h"

#include "Registration/cmtkImagePairSimilarityMeasureCR.h"
#include "Registration/cmtkImagePairSimilarityMeasureMSD.h"
#include "Registration/cmtkImagePairSimilarityMeasureNCC.h"
#include "Registration/cmtkImagePairSimilarityMeasureNMI.h"
#include "Registration/cmtkImagePairSimilarityMeasureMI.h"

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairNonrigidRegistrationFunctional::ImagePairNonrigidRegistrationFunctional
( UniformVolume::SmartPtr& reference, UniformVolume::SmartPtr& floating )
  : ImagePairRegistrationFunctional( reference, floating )
{
  this->m_NumberOfThreads = ThreadPool::GetGlobalThreadPool().GetNumberOfThreads();
  this->m_NumberOfTasks = 4 * this->m_NumberOfThreads - 3;
    
  Dim = 0;

  ReferenceFrom = UniformVolume::CoordinateVectorType( UniformVolume::CoordinateVectorType::Init( 0 ) );
  ReferenceTo = reference->Size;

  this->m_AdaptiveFixParameters = false;
  this->m_AdaptiveFixThreshFactor = 0.5;
  this->VolumeOfInfluence = NULL;

  this->m_ThreadWarp = Memory::AllocateArray<SplineWarpXform::SmartPtr>( this->m_NumberOfThreads );  

  this->m_ThreadVectorCache = Memory::AllocateArray<Vector3D*>( this->m_NumberOfThreads );
  for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
    this->m_ThreadVectorCache[thread] = Memory::AllocateArray<Vector3D>( this->m_ReferenceDims[0] );

  this->m_WarpedVolume = NULL;
  
  this->m_DimsX = this->m_ReferenceGrid->GetDims()[0];
  this->m_DimsY = this->m_ReferenceGrid->GetDims()[1];
  this->m_DimsZ = this->m_ReferenceGrid->GetDims()[2];
  
  this->m_FltDimsX = this->m_FloatingGrid->GetDims()[0];
  this->m_FltDimsY = this->m_FloatingGrid->GetDims()[1];
}

ImagePairNonrigidRegistrationFunctional::~ImagePairNonrigidRegistrationFunctional()
{
  for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread )
    if ( this->m_ThreadVectorCache[thread] ) 
      Memory::DeleteArray( this->m_ThreadVectorCache[thread] );
  Memory::DeleteArray( this->m_ThreadVectorCache );
  
  Memory::DeleteArray( this->m_ThreadWarp );
}

void
ImagePairNonrigidRegistrationFunctional::WeightedDerivative
( double& lower, double& upper, SplineWarpXform& warp, 
  const int param, const Types::Coordinate step ) const
{
  if ( this->m_JacobianConstraintWeight > 0 )
    {
    double lowerConstraint = 0, upperConstraint = 0;
    warp.GetJacobianConstraintDerivative( lowerConstraint, upperConstraint, param, VolumeOfInfluence[param], step );
    lower -= this->m_JacobianConstraintWeight * lowerConstraint;
    upper -= this->m_JacobianConstraintWeight * upperConstraint;
    } 

  if ( this->m_GridEnergyWeight > 0 ) 
    {
    double lowerEnergy = 0, upperEnergy = 0;
    warp.GetGridEnergyDerivative( lowerEnergy, upperEnergy, param, step );
    lower -= this->m_GridEnergyWeight * lowerEnergy;
    upper -= this->m_GridEnergyWeight * upperEnergy;
    }

  // Catch infinite values that result from a folding grid. Effectively
  // prevent this by setting the gradient term to 0.
  if ( !finite(upper) || !finite(lower) ) 
    {
    lower = upper = 0;
    }
  else
    {
    if ( this->m_MatchedLandmarkList ) 
      {
      double lowerMSD, upperMSD;
      warp.GetDerivativeLandmarksMSD( lowerMSD, upperMSD, this->m_MatchedLandmarkList, param, step );
      lower -= this->m_LandmarkErrorWeight * lowerMSD;
      upper -= this->m_LandmarkErrorWeight * upperMSD;
      }
    if ( this->m_InverseTransformation ) 
      {
      double lowerIC, upperIC;
      warp.GetDerivativeInverseConsistencyError( lowerIC, upperIC, this->m_InverseTransformation, this->m_ReferenceGrid, &(this->VolumeOfInfluence[param]), param, step );
      lower -= this->m_InverseConsistencyWeight * lowerIC;
      upper -= this->m_InverseConsistencyWeight * upperIC;
      }
    }
}

void
ImagePairNonrigidRegistrationFunctional::SetWarpXform
( SplineWarpXform::SmartPtr& warp )
{
  this->m_Warp = warp;
  if ( this->m_Warp )
    {
    this->m_Warp->RegisterVolume( this->m_ReferenceGrid );
    if ( Dim != this->m_Warp->VariableParamVectorDim() ) 
      {
      Dim = this->m_Warp->VariableParamVectorDim();
      this->m_StepScaleVector.resize( Dim );
      VolumeOfInfluence = Memory::AllocateArray<DataGrid::RegionType>( Dim );
      }
    
    DataGrid::RegionType *VOIptr = VolumeOfInfluence;
    Vector3D fromVOI, toVOI;
    for ( size_t dim=0; dim<Dim; ++dim, ++VOIptr ) 
      {
      this->m_StepScaleVector[dim] = this->GetParamStep( dim );
      this->m_Warp->GetVolumeOfInfluence( dim, ReferenceFrom, ReferenceTo, fromVOI, toVOI );
      *VOIptr = this->GetReferenceGridRange( fromVOI, toVOI );
      }
    
    for ( size_t thread = 0; thread < this->m_NumberOfThreads; ++thread ) 
      {
      if ( thread ) 
	{
	this->m_ThreadWarp[thread] = this->m_Warp->Clone();
	this->m_ThreadWarp[thread]->RegisterVolume( this->m_ReferenceGrid );
	} 
      else 
	{
	this->m_ThreadWarp[thread] = this->m_Warp;
	}
      } 
    }
}

ImagePairNonrigidRegistrationFunctional* 
ImagePairNonrigidRegistrationFunctional::Create
( const int metric,
  UniformVolume::SmartPtr& refVolume, 
  UniformVolume::SmartPtr& fltVolume,
  const Interpolators::InterpolationEnum interpolation )
{
  switch ( metric ) 
    {
    case 0:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNMI>( refVolume, fltVolume, interpolation );
    case 1:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMI>( refVolume, fltVolume, interpolation );
    case 2:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureCR>( refVolume, fltVolume, interpolation );
    case 3:
      return NULL; // masked NMI retired
    case 4:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureMSD>( refVolume, fltVolume, interpolation );
    case 5:
      return new ImagePairNonrigidRegistrationFunctionalTemplate<ImagePairSimilarityMeasureNCC>( refVolume, fltVolume, interpolation );
    default:
      return NULL;
    }

  return NULL;
}

} // namespace cmtk
