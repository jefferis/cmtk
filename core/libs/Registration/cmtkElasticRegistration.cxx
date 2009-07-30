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

#include <cmtkElasticRegistration.h>

#include <cmtkLandmarkList.h>
#include <cmtkMatchedLandmarkList.h>

#include <cmtkVoxelMatchingElasticFunctional.h>
#include <cmtkSymmetricElasticFunctional.h>

#include <cmtkOptimizer.h>
#include <cmtkBestNeighbourOptimizer.h>
#include <cmtkBestDirectionOptimizer.h>

#include <cmtkUniformVolume.h>
#include <cmtkSplineWarpXform.h>

#include <cmtkReformatVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ElasticRegistration::ElasticRegistration () 
  : VoxelRegistration(),
    InitialWarpXform( NULL ),
    InverseWarpXform( NULL ),
    m_RigidityConstraintMap( NULL ),
    m_InverseConsistencyWeight( 0.0 ),
    m_ForceOutsideFlag( false ),
    m_ForceOutsideValue( 0.0 )
{
  ForceSwitchVolumes = 0;
  this->m_GridSpacing = 10;
  RestrictToAxes = NULL;
  this->m_RefineGrid = 0;
  RefinedGridAtLevel = -1;
  RefineGridCount = 0;
  this->m_DelayRefineGrid = 0;
  RefineDelayed = false;
  IgnoreEdge = 0;
  this->m_FastMode = 0;
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
  CallbackResult result = this->Superclass::InitRegistration();
  if ( result != CALLBACK_OK ) return result;

  UniformVolume::SmartPtr refVolume;
  UniformVolume::SmartPtr fltVolume;

  if ( ForceSwitchVolumes ) 
    {
    refVolume = this->m_Volume_2;
    fltVolume = this->m_Volume_1;
    SwitchVolumes = 1;
    } 
  else
    {
    refVolume = this->m_Volume_1;
    fltVolume = this->m_Volume_2;
    SwitchVolumes = 0;
    }

  MatchedLandmarkList::SmartPtr mll( NULL );
  if ( this->m_LandmarkErrorWeight != 0 ) 
    {
    LandmarkList::SmartPtr sourceLandmarks = refVolume->m_LandmarkList;
    LandmarkList::SmartPtr targetLandmarks = fltVolume->m_LandmarkList;
    
    if ( sourceLandmarks && targetLandmarks ) 
      {
      mll = MatchedLandmarkList::SmartPtr( new MatchedLandmarkList( sourceLandmarks, targetLandmarks ) );
      fprintf( stderr, "Matched %d landmarks.\n", (int)mll->size() );
      }
    }
  
  AffineXform::SmartPtr affineXform;
  AffineXform::SmartPtr initialInverse;
  if ( SwitchVolumes ) 
    {
    initialInverse = AffineXform::SmartPtr( static_cast<AffineXform*>( this->m_InitialXform->MakeInverse() ) );
    CoordinateVector v( initialInverse->ParamVectorDim() );
    affineXform = AffineXform::SmartPtr( new AffineXform( initialInverse->GetParamVector(v) ) );
    initialInverse = this->m_InitialXform;
    } 
  else
    {
    affineXform = this->m_InitialXform;
    initialInverse = AffineXform::SmartPtr( static_cast<AffineXform*>( this->m_InitialXform->MakeInverse() ) );
    }
  
  Vector3D center = fltVolume->GetCenterCropRegion();
  affineXform->ChangeCenter( center.XYZ );

  Types::Coordinate currSampling = std::max( this->m_Sampling, 2 * std::min( refVolume->GetMinDelta(), fltVolume->GetMinDelta()));

  // If no initial transformation exists, create one from the defined
  // parameters.
  if ( InitialWarpXform ) 
    {
    // If we have an initial transformation from somewhere, use that.
    // This will override all registration parameters otherwise defined,
    // for example grid spacing and deformation type.
    InitialWarpXform->SetIgnoreEdge( IgnoreEdge );
    InitialWarpXform->SetFastMode( this->m_FastMode );
    InitialWarpXform->SetParametersActive( RestrictToAxes );
    // MIPSpro needs explicit.
    this->m_Xform = Xform::SmartPtr::DynamicCastFrom( InitialWarpXform );
    } 
  else
    {
    WarpXform::SmartPtr warpXform( this->MakeWarpXform( refVolume->Size, initialInverse ) );
    
    if ( this->m_InverseConsistencyWeight > 0 ) 
      InverseWarpXform = WarpXform::SmartPtr( this->MakeWarpXform( fltVolume->Size, affineXform ) );

    // MIPSpro needs explicit:
    this->m_Xform = Xform::SmartPtr::DynamicCastFrom( warpXform ); 
    }
  
  if ( this->m_UseOriginalData )
    {
    Functional::SmartPtr nextFunctional( this->MakeFunctional( refVolume, fltVolume, this->m_RigidityConstraintMap, mll ) );
    FunctionalStack.push( nextFunctional );
    }
  
  CallbackResult irq = CALLBACK_OK;
  
  double coarsest = CoarsestResolution;
  if ( coarsest <= 0 ) 
    coarsest = this->m_Exploration;
  
  for ( ;(currSampling<=coarsest) && ( irq == CALLBACK_OK ); currSampling *= 2 ) 
    {
    UniformVolume::SmartPtr nextRef;
    UniformVolume::SmartPtr nextMod;
    UniformVolume::SmartPtr nextRigidityMap;
    try 
      {
      this->ReportProgress( "resampling", 0 );
      nextRef = UniformVolume::SmartPtr( new UniformVolume( *refVolume, currSampling ) );
      irq = this->ReportProgress(  "resampling", 50 );
      nextMod = UniformVolume::SmartPtr( new UniformVolume( *fltVolume, currSampling ) );
      if ( this->m_RigidityConstraintMap )
	{
	nextRigidityMap = UniformVolume::SmartPtr( new UniformVolume( *this->m_RigidityConstraintMap, currSampling ) );
	}
      }
    catch (...) 
      {
      }
    
    Functional::SmartPtr nextFunctional( this->MakeFunctional( nextRef, nextMod, nextRigidityMap, mll ) );
    FunctionalStack.push( nextFunctional );
    
    refVolume = nextRef;
    fltVolume = nextMod;
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
  return irq;
}

WarpXform*
ElasticRegistration::MakeWarpXform
( const Types::Coordinate* size, const AffineXform* initialAffine ) const
{
  WarpXform* warpXform = NULL;
  
  warpXform = new SplineWarpXform( size, this->m_GridSpacing, initialAffine, this->m_ExactGridSpacing );
  
  warpXform->SetIgnoreEdge( this->IgnoreEdge );
  warpXform->SetFastMode( this->m_FastMode );
  warpXform->SetParametersActive( this->RestrictToAxes );

  return warpXform;
}

Functional* 
ElasticRegistration::MakeFunctional
( UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume,
  UniformVolume::SmartPtr& rigidityMap,
  MatchedLandmarkList::SmartPtr& mll ) const
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
    if ( rigidityMap )
      {
      newFunctional->SetRigidityConstraintMap( rigidityMap );
      }
    newFunctional->SetGridEnergyWeight( this->m_GridEnergyWeight );
    if ( mll )
      {
      newFunctional->SetLandmarkErrorWeight( this->m_LandmarkErrorWeight );
      newFunctional->SetMatchedLandmarkList( mll );
      }
    
    return newFunctional;
  }
}

void
ElasticRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& functional,
  const int idx, const int total ) 
{
  WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->m_Xform );

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
	  warpXform->Refine( 2 );
	  if ( InverseWarpXform )
	    InverseWarpXform->Refine( 2 );
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

UniformVolume* 
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
