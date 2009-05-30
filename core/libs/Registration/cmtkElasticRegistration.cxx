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

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ElasticRegistration::ElasticRegistration () 
  : VoxelRegistration(),
    InitialWarpXform( NULL ),
    InverseWarpXform( NULL ),
    RigidityConstraintMap( NULL ),
    InverseConsistencyWeight( 0.0 ),
    m_ForceOutsideFlag( false ),
    m_ForceOutsideValue( 0.0 )
{
  ForceSwitchVolumes = 0;
  GridSpacing = 10;
  RestrictToAxes = NULL;
  RefineGrid = 0;
  RefinedGridAtLevel = -1;
  RefineGridCount = 0;
  DelayRefineGrid = 0;
  RefineDelayed = false;
  IgnoreEdge = 0;
  FastMode = 0;
  AdaptiveFixParameters = 1;
  AdaptiveFixThreshFactor = 0.5;
  JacobianConstraintWeight = 0;
  RigidityConstraintWeight = 0;
  GridEnergyWeight = 0;
  RelaxWeight = -1;
  LandmarkErrorWeight = 0;
  InverseConsistencyWeight = 0.0;
  RelaxationStep = false;
  Sobel1 = Sobel2 = 0;
  HistogramEqualization1 = HistogramEqualization2 = 0;
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
    refVolume = this->Volume_2;
    fltVolume = this->Volume_1;
    SwitchVolumes = 1;
    } 
  else
    {
    refVolume = this->Volume_1;
    fltVolume = this->Volume_2;
    SwitchVolumes = 0;
    }

  MatchedLandmarkList::SmartPtr mll( NULL );
  if ( LandmarkErrorWeight != 0 ) 
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
    initialInverse = AffineXform::SmartPtr( static_cast<AffineXform*>( InitialXform->MakeInverse() ) );
    CoordinateVector v( initialInverse->ParamVectorDim() );
    affineXform = AffineXform::SmartPtr( new AffineXform( initialInverse->GetParamVector(v) ) );
    initialInverse = InitialXform;
    } 
  else
    {
    affineXform = InitialXform;
    initialInverse = AffineXform::SmartPtr( static_cast<AffineXform*>( InitialXform->MakeInverse() ) );
    }
  
  Vector3D center = fltVolume->GetCenterCropRegion();
  affineXform->ChangeCenter( center.XYZ );

  Types::Coordinate currSampling = std::max( Sampling, 2 * std::min( refVolume->GetMinDelta(), fltVolume->GetMinDelta()));

  // If no initial transformation exists, create one from the defined
  // parameters.
  if ( InitialWarpXform ) 
    {
    // If we have an initial transformation from somewhere, use that.
    // This will override all registration parameters otherwise defined,
    // for example grid spacing and deformation type.
    InitialWarpXform->SetIgnoreEdge( IgnoreEdge );
    InitialWarpXform->SetFastMode( FastMode );
    InitialWarpXform->SetParametersActive( RestrictToAxes );
    // MIPSpro needs explicit.
    Xform = Xform::SmartPtr::DynamicCastFrom( InitialWarpXform );
    } 
  else
    {
    WarpXform::SmartPtr warpXform( this->MakeWarpXform( refVolume->Size, initialInverse ) );
    
    if ( this->InverseConsistencyWeight > 0 ) 
      InverseWarpXform = WarpXform::SmartPtr( this->MakeWarpXform( fltVolume->Size, affineXform ) );

    // MIPSpro needs explicit:
    Xform = Xform::SmartPtr::DynamicCastFrom( warpXform ); 
    }
  
  if ( UseOriginalData )
    {
    Functional::SmartPtr nextFunctional( this->MakeFunctional( refVolume, fltVolume, this->RigidityConstraintMap, mll ) );
    FunctionalStack.push( nextFunctional );
    }
  
  CallbackResult irq = CALLBACK_OK;
  
  double coarsest = CoarsestResolution;
  if ( coarsest <= 0 ) coarsest = Exploration;
  
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
      if ( this->RigidityConstraintMap )
	{
	nextRigidityMap = UniformVolume::SmartPtr( new UniformVolume( *this->RigidityConstraintMap, currSampling ) );
	}
      }
    catch (...) 
      {
      }
    
    UniformVolume::SmartPtr useRef = nextRef;
    UniformVolume::SmartPtr useMod = nextMod;
    UniformVolume::SmartPtr useRigidityMap = nextRigidityMap;
    
    if ( HistogramEqualization1 ) 
      {
      useRef = UniformVolume::SmartPtr( useRef->Clone() );
      useRef->GetData()->HistogramEqualization();
      }
    if ( HistogramEqualization2 )
      {
      useMod = UniformVolume::SmartPtr( useMod->Clone() );
      useMod->GetData()->HistogramEqualization();
      }

    if ( Sobel1 ) 
      {
      useRef = UniformVolume::SmartPtr( useRef->Clone() );
      useRef->ApplySobelFilter();
      }
    if ( Sobel2 ) 
      {
      useMod = UniformVolume::SmartPtr( useMod->Clone() );
      useMod->ApplySobelFilter();
      }

    Functional::SmartPtr nextFunctional
      ( this->MakeFunctional( useRef, useMod, useRigidityMap, mll ) );
    FunctionalStack.push( nextFunctional );
    
    refVolume = nextRef;
    fltVolume = nextMod;
  }
  
  switch ( Algorithm ) {
  case 0:
    Optimizer = Optimizer::SmartPtr
      ( new BestNeighbourOptimizer( OptimizerStepFactor ) ); 
    break;
  case 1:
  case 2:
    Optimizer = Optimizer::SmartPtr( NULL );
    break;
  case 3: {
    BestDirectionOptimizer *optimizer = 
      new BestDirectionOptimizer( OptimizerStepFactor ); 
    optimizer->SetUseMaxNorm( UseMaxNorm );
    Optimizer = Optimizer::SmartPtr( optimizer );
    break;
  }
  }

  Optimizer->SetCallback( Callback );
  return irq;
}

WarpXform*
ElasticRegistration::MakeWarpXform
( const Types::Coordinate* size, const AffineXform* initialAffine ) const
{
  WarpXform* warpXform = NULL;
  
  warpXform = new SplineWarpXform( size, this->GridSpacing, initialAffine, this->ExactGridSpacing );
  
  warpXform->SetIgnoreEdge( this->IgnoreEdge );
  warpXform->SetFastMode( this->FastMode );
  warpXform->SetParametersActive( this->RestrictToAxes );

  return warpXform;
}

Functional* 
ElasticRegistration::MakeFunctional
( UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume,
  UniformVolume::SmartPtr& rigidityMap,
  MatchedLandmarkList::SmartPtr& mll ) const
{
  if ( this->InverseConsistencyWeight > 0 ) 
    {
    SymmetricElasticFunctional *newFunctional = CreateSymmetricElasticFunctional( this->Metric, refVolume, fltVolume );
    newFunctional->SetInverseConsistencyWeight( this->InverseConsistencyWeight );
    newFunctional->SetAdaptiveFixParameters( this->AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->JacobianConstraintWeight );
    newFunctional->SetRigidityConstraintWeight( this->RigidityConstraintWeight );
    newFunctional->SetGridEnergyWeight( this->GridEnergyWeight );
//    newFunctional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    return newFunctional;
    } 
  else
    {
    VoxelMatchingElasticFunctional *newFunctional = CreateElasticFunctional( this->Metric, refVolume, fltVolume );
    newFunctional->SetAdaptiveFixParameters( this->AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->JacobianConstraintWeight );
    newFunctional->SetRigidityConstraintWeight( this->RigidityConstraintWeight );
    newFunctional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    if ( rigidityMap )
      {
      newFunctional->SetRigidityConstraintMap( rigidityMap );
      }
    newFunctional->SetGridEnergyWeight( this->GridEnergyWeight );    
    if ( mll )
      {
      newFunctional->SetLandmarkErrorWeight( this->LandmarkErrorWeight );
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
  WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->Xform );

  float effGridEnergyWeight = this->GridEnergyWeight;
  float effJacobianConstraintWeight = this->JacobianConstraintWeight;
  float effRigidityConstraintWeight = this->RigidityConstraintWeight;
  float effInverseConsistencyWeight = this->InverseConsistencyWeight;

  if ( (this->RelaxWeight > 0) && !this->RelaxationStep ) {
    effGridEnergyWeight *= this->RelaxWeight;
    effJacobianConstraintWeight *= this->RelaxWeight;
    effRigidityConstraintWeight *= this->RelaxWeight;
    effInverseConsistencyWeight *= this->RelaxWeight;
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
  if ( ( RelaxWeight > 0 ) && ! RelaxationStep ) 
    {
    RelaxationStep = true;
    this->Superclass::DoneResolution( v, functional, idx, total );
    return false; // repeat with a relaxation step.
    } 
  else 
    {
    RelaxationStep = false;
    }
  
  bool repeat = ( ( idx == total ) && ( RefineGridCount < RefineGrid ) );
  
  if ( (RefinedGridAtLevel != idx) || (idx==total) ) 
    {    
    if ( RefineGridCount < RefineGrid ) 
      {      
      if ( (!DelayRefineGrid) || RefineDelayed || ( idx == total ) ) 
	{
	WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( Xform );
	if ( warpXform ) 
	  {
	  warpXform->Refine( 2 );
	  if ( InverseWarpXform )
	    InverseWarpXform->Refine( 2 );
	  ++RefineGridCount;
	  functional->GetParamVector( *v );    
	  if ( Callback ) 
	    Callback->Comment( "Refined control point grid." );
	  RefinedGridAtLevel = idx;
	  } 	  
	if ( DelayRefineGrid && ( idx > 1 ) ) repeat = true;
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

} // namespace cmtk
