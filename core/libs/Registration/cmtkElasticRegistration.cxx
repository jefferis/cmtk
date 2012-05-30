/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include "cmtkElasticRegistration.h"

#include <Base/cmtkLandmarkList.h>
#include <Base/cmtkLandmarkPairList.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>

#include <Registration/cmtkVoxelMatchingElasticFunctional.h>
#include <Registration/cmtkSymmetricElasticFunctional.h>
#include <Registration/cmtkOptimizer.h>
#include <Registration/cmtkBestNeighbourOptimizer.h>
#include <Registration/cmtkBestDirectionOptimizer.h>
#include <Registration/cmtkReformatVolume.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ElasticRegistration::ElasticRegistration () 
  : VoxelRegistration(),
    InitialWarpXform( NULL ),
    InverseWarpXform( NULL ),
    ForceSwitchVolumes( false ),
    m_MatchFltToRefHistogram( false ),
    m_RigidityConstraintMap( NULL ),
    m_InverseConsistencyWeight( 0.0 ),
    m_RelaxToUnfold( false ),
    m_ForceOutsideFlag( false ),
    m_ForceOutsideValue( 0.0 )
{
  this->m_GridSpacing = 10;
  RestrictToAxes = NULL;
  this->m_RefineGrid = 0;
  RefinedGridAtLevel = -1;
  RefineGridCount = 0;
  this->m_DelayRefineGrid = 0;
  RefineDelayed = false;
  IgnoreEdge = 0;
  this->m_FastMode = false;
  this->m_AdaptiveFixParameters = 1;
  this->m_AdaptiveFixThreshFactor = 0.5;
  this->m_JacobianConstraintWeight = 0;
  this->m_RigidityConstraintWeight = 0;
  this->m_GridEnergyWeight = 0;
  this->m_RelaxWeight = -1;
  this->m_LandmarkErrorWeight = 0;
  this->m_InverseConsistencyWeight = 0.0;
  RelaxationStep = false;
}

CallbackResult 
ElasticRegistration::InitRegistration ()
{
  this->m_ReferenceVolume = this->m_Volume_1;
  this->m_FloatingVolume = this->m_Volume_2;

  if ( this->m_MatchFltToRefHistogram )
    {
    this->GetVolume_2()->GetData()->ApplyFunctionObject( TypedArrayFunctionHistogramMatching( *(this->GetVolume_2()->GetData()), *(this->GetVolume_1()->GetData()) ) );
    }
  
  if ( this->m_LandmarkErrorWeight != 0 ) 
    {
    LandmarkList::SmartPtr sourceLandmarks = this->m_ReferenceVolume->m_LandmarkList;
    LandmarkList::SmartPtr targetLandmarks = this->m_FloatingVolume->m_LandmarkList;
    
    if ( sourceLandmarks && targetLandmarks ) 
      {
      this->m_LandmarkPairs = LandmarkPairList::SmartPtr( new LandmarkPairList( *(sourceLandmarks), *(targetLandmarks) ) );
      StdErr << "Matched " << this->m_LandmarkPairs->size() << " landmarks.\n";
      }
    }
  
  AffineXform::SmartPtr affineXform = this->m_InitialTransformation;
  AffineXform::SmartPtr initialInverse = AffineXform::SmartPtr::DynamicCastFrom( this->m_InitialTransformation->GetInverse() );
  
  affineXform->ChangeCenter( this->m_FloatingVolume->GetCenterCropRegion() );

  Types::Coordinate currSampling = std::max( this->m_Sampling, 2 * std::min( this->m_ReferenceVolume->GetMinDelta(), this->m_FloatingVolume->GetMinDelta()));

  // If no initial transformation exists, create one from the defined
  // parameters.
  if ( InitialWarpXform ) 
    {
    // If we have an initial transformation from somewhere, use that.
    // This will override all registration parameters otherwise defined,
    // for example grid spacing and deformation type.
    InitialWarpXform->SetIgnoreEdge( IgnoreEdge );
    InitialWarpXform->SetFastMode( this->m_FastMode );
    // MIPSpro needs explicit.
    this->m_Xform = Xform::SmartPtr::DynamicCastFrom( InitialWarpXform );
    } 
  else
    {
    SplineWarpXform::SmartPtr warpXform( this->MakeWarpXform( this->m_ReferenceVolume->m_Size, affineXform ) );
    
    if ( this->m_InverseConsistencyWeight > 0 ) 
      InverseWarpXform = SplineWarpXform::SmartPtr( this->MakeWarpXform( this->m_FloatingVolume->m_Size, initialInverse ) );

    // MIPSpro needs explicit:
    this->m_Xform = Xform::SmartPtr::DynamicCastFrom( warpXform ); 
    }
  
  if ( this->m_UseOriginalData )
    {
    Functional::SmartPtr nextFunctional( this->MakeFunctional( this->m_ReferenceVolume, this->m_FloatingVolume, this->m_RigidityConstraintMap ) );
    FunctionalStack.push( nextFunctional );
    }
  
  if ( this->m_Exploration <= 0 )
    {
    const SplineWarpXform* warp = SplineWarpXform::SmartPtr::DynamicCastFrom( this->m_Xform ); 
    this->m_Exploration = 0.25 * std::max( warp->m_Spacing[0], std::max( warp->m_Spacing[1], warp->m_Spacing[2] ) );
    }

  if ( this->CoarsestResolution <= 0 ) 
    this->CoarsestResolution = this->m_Exploration;
  
  UniformVolume::SmartPtr currRef( this->m_ReferenceVolume );
  UniformVolume::SmartPtr currFlt( this->m_FloatingVolume );

  for ( ;(currSampling<=this->CoarsestResolution); currSampling *= 2 ) 
    {
    UniformVolume::SmartPtr nextRef( new UniformVolume( *currRef, currSampling ) );
    UniformVolume::SmartPtr nextFlt( new UniformVolume( *currFlt, currSampling ) );

    UniformVolume::SmartPtr nextRigidityMap;
    if ( this->m_RigidityConstraintMap )
      {
      nextRigidityMap = UniformVolume::SmartPtr( new UniformVolume( *this->m_RigidityConstraintMap, currSampling ) );
      }
    
    Functional::SmartPtr nextFunctional( this->MakeFunctional( nextRef, nextFlt, nextRigidityMap ) );
    FunctionalStack.push( nextFunctional );
    
    currRef = nextRef;
    currFlt = nextFlt;
    }
  
  switch ( this->m_Algorithm ) 
    {
    case 0:
      this->m_Optimizer = Optimizer::SmartPtr( new BestNeighbourOptimizer( OptimizerStepFactor ) ); 
      break;
    case 1:
    case 2:
      this->m_Optimizer = Optimizer::SmartPtr( NULL );
      break;
    case 3: 
    {
    BestDirectionOptimizer *optimizer = new BestDirectionOptimizer( OptimizerStepFactor ); 
    optimizer->SetUseMaxNorm( UseMaxNorm );
    this->m_Optimizer = Optimizer::SmartPtr( optimizer );
    break;
    }
    }
  
  this->m_Optimizer->SetCallback( this->m_Callback );

  return this->Superclass::InitRegistration();
}

const SplineWarpXform::SmartPtr
ElasticRegistration::MakeWarpXform
( const Vector3D& size, const AffineXform* initialAffine ) const
{
  SplineWarpXform::SmartPtr warpXform( new SplineWarpXform( size, this->m_GridSpacing, initialAffine, this->m_ExactGridSpacing ) );
  
  warpXform->SetIgnoreEdge( this->IgnoreEdge );
  warpXform->SetFastMode( this->m_FastMode );
  warpXform->SetParametersActive();

  return warpXform;
}

Functional* 
ElasticRegistration::MakeFunctional
( UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume,
  UniformVolume::SmartPtr& rigidityMap ) const
{
  if ( this->m_InverseConsistencyWeight > 0 ) 
    {
    SymmetricElasticFunctional *newFunctional = CreateSymmetricElasticFunctional( this->m_Metric, refVolume, fltVolume );
    newFunctional->SetInverseConsistencyWeight( this->m_InverseConsistencyWeight );
    newFunctional->SetAdaptiveFixParameters( this->m_AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->m_AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->m_JacobianConstraintWeight );
    newFunctional->SetRigidityConstraintWeight( this->m_RigidityConstraintWeight );
    newFunctional->SetGridEnergyWeight( this->m_GridEnergyWeight );
//    newFunctional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    return newFunctional;
    } 
  else
    {
    VoxelMatchingElasticFunctional *newFunctional = CreateElasticFunctional( this->m_Metric, refVolume, fltVolume );
    newFunctional->SetAdaptiveFixParameters( this->m_AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->m_AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->m_JacobianConstraintWeight );
    newFunctional->SetRigidityConstraintWeight( this->m_RigidityConstraintWeight );
    newFunctional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    newFunctional->SetActiveCoordinates( this->RestrictToAxes );
    if ( rigidityMap )
      {
      newFunctional->SetRigidityConstraintMap( rigidityMap );
      }
    newFunctional->SetGridEnergyWeight( this->m_GridEnergyWeight );
    if ( this->m_LandmarkPairs )
      {
      newFunctional->SetLandmarkErrorWeight( this->m_LandmarkErrorWeight );
      newFunctional->SetLandmarkPairs( this->m_LandmarkPairs );
      }
    
    return newFunctional;
  }
}

void
ElasticRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& functional,
  const int idx, const int total ) 
{
  SplineWarpXform::SmartPtr warpXform = SplineWarpXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  
  float effGridEnergyWeight = this->m_GridEnergyWeight;
  float effJacobianConstraintWeight = this->m_JacobianConstraintWeight;
  float effRigidityConstraintWeight = this->m_RigidityConstraintWeight;
  float effInverseConsistencyWeight = this->m_InverseConsistencyWeight;

  if ( (this->m_RelaxWeight > 0) && !this->RelaxationStep ) 
    {
    effGridEnergyWeight *= this->m_RelaxWeight;
    effJacobianConstraintWeight *= this->m_RelaxWeight;
    effRigidityConstraintWeight *= this->m_RelaxWeight;
    effInverseConsistencyWeight *= this->m_RelaxWeight;
    }

  // handle simple elastic functional
  SmartPointer<VoxelMatchingElasticFunctional> elasticFunctional = VoxelMatchingElasticFunctional::SmartPtr::DynamicCastFrom( functional );
  if ( elasticFunctional ) 
    {
    elasticFunctional->SetWarpXform( warpXform );

    if ( this->m_RelaxToUnfold )
      warpXform->RelaxToUnfold();

    elasticFunctional->SetGridEnergyWeight( effGridEnergyWeight );
    elasticFunctional->SetJacobianConstraintWeight( effJacobianConstraintWeight );
    elasticFunctional->SetRigidityConstraintWeight( effRigidityConstraintWeight );
    } 
  else 
    {
    // handle inverse-consistent elastic functional
    SmartPointer<SymmetricElasticFunctional> symmetricFunctional = SymmetricElasticFunctional::SmartPtr::DynamicCastFrom( functional );
    if ( symmetricFunctional ) 
      {
      symmetricFunctional->SetWarpXform( warpXform, this->InverseWarpXform );

      if ( this->m_RelaxToUnfold )
	{
	warpXform->RelaxToUnfold();
	this->InverseWarpXform->RelaxToUnfold();
	}

      symmetricFunctional->SetGridEnergyWeight( effGridEnergyWeight );
      symmetricFunctional->SetJacobianConstraintWeight( effJacobianConstraintWeight );
      symmetricFunctional->SetRigidityConstraintWeight( effRigidityConstraintWeight );
      symmetricFunctional->SetInverseConsistencyWeight( effInverseConsistencyWeight );
      } 
    else 
      {
      // neither simple nor inverse-consistent functional: something went
      // badly wrong.
      StdErr << "Fatal coding error: Non-elastic functional in ElasticRegistration::EnterResolution.\n";
      abort();
      }
    }
  
  Superclass::EnterResolution( v, functional, idx, total );
}

int 
ElasticRegistration::DoneResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& functional,
  const int idx, const int total ) 
{
  if ( ( this->m_RelaxWeight > 0 ) && ! RelaxationStep ) 
    {
    RelaxationStep = true;
    this->Superclass::DoneResolution( v, functional, idx, total );
    return false; // repeat with a relaxation step.
    } 
  else 
    {
    RelaxationStep = false;
    }
  
  bool repeat = ( ( idx == total ) && ( RefineGridCount < this->m_RefineGrid ) );
  
  if ( (RefinedGridAtLevel != idx) || (idx==total) ) 
    {    
    if ( RefineGridCount < this->m_RefineGrid ) 
      {      
      if ( (!this->m_DelayRefineGrid) || RefineDelayed || ( idx == total ) ) 
	{
	WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->m_Xform );
	if ( warpXform ) 
	  {
	  warpXform->Refine();
	  if ( InverseWarpXform )
	    InverseWarpXform->Refine();
	  ++RefineGridCount;
	  functional->GetParamVector( *v );    
	  if ( this->m_Callback ) 
	    this->m_Callback->Comment( "Refined control point grid." );
	  RefinedGridAtLevel = idx;
	  } 	  
	if ( this->m_DelayRefineGrid && ( idx > 1 ) ) repeat = true;
	RefineDelayed = false;
	} 
      else 
	{
	RefineDelayed = true;
	}
      }
    }
  else 
    {
    RefineDelayed = true;
    }
  
  return this->Superclass::DoneResolution( v, functional, idx, total ) && !repeat;
}

const UniformVolume::SmartPtr
ElasticRegistration::GetReformattedFloatingImage( Interpolators::InterpolationEnum interpolator )
{
  ReformatVolume reformat;
  reformat.SetInterpolation( interpolator );
  reformat.SetReferenceVolume( this->m_Volume_1 );
  reformat.SetFloatingVolume( this->m_Volume_2 );

  WarpXform::SmartPtr warpXform( this->GetTransformation() );
  reformat.SetWarpXform( warpXform );

  return reformat.PlainReformat();
}


} // namespace cmtk
