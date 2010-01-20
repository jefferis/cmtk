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

#include <cmtkImagePairNonrigidRegistration.h>

#include <cmtkLandmarkList.h>
#include <cmtkMatchedLandmarkList.h>

#include <cmtkImagePairNonrigidRegistrationFunctional.h>
#include <cmtkImagePairSymmetricNonrigidRegistrationFunctional.h>

#include <cmtkOptimizer.h>
#include <cmtkBestNeighbourOptimizer.h>
#include <cmtkBestDirectionOptimizer.h>

#include <cmtkUniformVolume.h>
#include <cmtkSplineWarpXform.h>

#include <cmtkReformatVolume.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairNonrigidRegistration::ImagePairNonrigidRegistration () 
  : InitialWarpXform( NULL ),
    InverseWarpXform( NULL ),
    m_MatchFltToRefHistogram( false ),
    m_InverseConsistencyWeight( 0.0 ),
    m_ForceOutsideFlag( false ),
    m_ForceOutsideValue( 0.0 )
{
  this->m_Metric = 0;
  this->m_Algorithm = 3;

  this->m_GridSpacing = 15;
  this->m_ExactGridSpacing = 0;
  this->m_GridSpacing = 10;
  RestrictToAxes = NULL;
  this->m_RefineGrid = 0;
  RefinedGridAtLevel = -1;
  RefineGridCount = 0;
  this->m_DelayRefineGrid = 0;
  RefineDelayed = false;
  IgnoreEdge = 0;
  this->m_FastMode = true;
  this->m_AdaptiveFixParameters = 1;
  this->m_AdaptiveFixThreshFactor = 0.5;
  this->m_JacobianConstraintWeight = 0;
  this->m_GridEnergyWeight = 0;
  this->m_RelaxWeight = -1;
  this->m_LandmarkErrorWeight = 0;
  this->m_InverseConsistencyWeight = 0.0;
  RelaxationStep = false;
}

CallbackResult 
ImagePairNonrigidRegistration::InitRegistration ()
{
  this->m_ReferenceVolume = this->m_Volume_1;
  this->m_FloatingVolume = this->m_Volume_2;

  if ( this->m_MatchFltToRefHistogram )
    {
    this->GetVolume_2()->GetData()->MatchHistogramToReference( this->GetVolume_1()->GetData() );
    }
  
  MatchedLandmarkList::SmartPtr mll( NULL );
  if ( this->m_LandmarkErrorWeight != 0 ) 
    {
    LandmarkList::SmartPtr sourceLandmarks = this->m_ReferenceVolume->m_LandmarkList;
    LandmarkList::SmartPtr targetLandmarks = this->m_FloatingVolume->m_LandmarkList;
    
    if ( sourceLandmarks && targetLandmarks ) 
      {
      mll = MatchedLandmarkList::SmartPtr( new MatchedLandmarkList( sourceLandmarks, targetLandmarks ) );
      fprintf( stderr, "Matched %d landmarks.\n", (int)mll->size() );
      }
    }
  
  Vector3D center = this->m_FloatingVolume->GetCenterCropRegion();
  this->m_InitialTransformation->ChangeCenter( center.XYZ );

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
    InitialWarpXform->SetParametersActive( RestrictToAxes );
    // MIPSpro needs explicit.
    this->m_Xform = Xform::SmartPtr::DynamicCastFrom( InitialWarpXform );
    } 
  else
    {
    WarpXform::SmartPtr warpXform( this->MakeWarpXform( this->m_ReferenceVolume->Size, this->m_InitialTransformation ) );
    
    if ( this->m_InverseConsistencyWeight > 0 ) 
      InverseWarpXform = WarpXform::SmartPtr( this->MakeWarpXform( this->m_FloatingVolume->Size, this->m_InitialTransformation->GetInverse() ) );

    // MIPSpro needs explicit:
    this->m_Xform = Xform::SmartPtr::DynamicCastFrom( warpXform ); 
    }
  
  if ( this->m_UseOriginalData )
    {
    Functional::SmartPtr nextFunctional( this->MakeFunctional( this->m_ReferenceVolume, this->m_FloatingVolume, mll ) );
    FunctionalStack.push( nextFunctional );
    }
  
  if ( this->m_MaxStepSize <= 0 )
    {
    const SplineWarpXform* warp = SplineWarpXform::SmartPtr::DynamicCastFrom( this->m_Xform ); 
    this->m_MaxStepSize = 0.25 * std::max( warp->Spacing[0], std::max( warp->Spacing[1], warp->Spacing[2] ) );
    }

  if ( this->m_CoarsestResolution <= 0 ) 
    this->m_CoarsestResolution = this->m_MaxStepSize;
  
  UniformVolume::SmartPtr currRef( this->m_ReferenceVolume );
  UniformVolume::SmartPtr currFlt( this->m_FloatingVolume );

  for ( ;(currSampling<=this->m_CoarsestResolution); currSampling *= 2 ) 
    {
    UniformVolume::SmartPtr nextRef( new UniformVolume( *currRef, currSampling ) );
    UniformVolume::SmartPtr nextFlt( new UniformVolume( *currFlt, currSampling ) );

    Functional::SmartPtr nextFunctional( this->MakeFunctional( nextRef, nextFlt, mll ) );
    FunctionalStack.push( nextFunctional );
    
    currRef = nextRef;
    currFlt = nextFlt;
    }
  
  switch ( this->m_Algorithm ) 
    {
    case 0:
      this->m_Optimizer = Optimizer::SmartPtr( new BestNeighbourOptimizer( this->m_OptimizerStepFactor ) ); 
      break;
    case 1:
    case 2:
      this->m_Optimizer = Optimizer::SmartPtr( NULL );
      break;
    case 3: 
    {
    BestDirectionOptimizer *optimizer = new BestDirectionOptimizer( this->m_OptimizerStepFactor );
    optimizer->SetUseMaxNorm( this->m_UseMaxNorm );
    this->m_Optimizer = Optimizer::SmartPtr( optimizer );
    break;
    }
    }
  
  this->m_Optimizer->SetCallback( this->m_Callback );

  return this->Superclass::InitRegistration();
}

WarpXform*
ImagePairNonrigidRegistration::MakeWarpXform
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
ImagePairNonrigidRegistration::MakeFunctional
( UniformVolume::SmartPtr& refVolume, UniformVolume::SmartPtr& fltVolume,
  MatchedLandmarkList::SmartPtr& mll ) const
{
  WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  if ( this->m_InverseConsistencyWeight > 0 ) 
    {
    ImagePairSymmetricNonrigidRegistrationFunctional *newFunctional = 
      ImagePairSymmetricNonrigidRegistrationFunctional::Create( this->m_Metric, refVolume, fltVolume, this->m_FloatingImageInterpolation );
    newFunctional->SetInverseConsistencyWeight( this->m_InverseConsistencyWeight );
    newFunctional->SetAdaptiveFixParameters( this->m_AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->m_AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->m_JacobianConstraintWeight );
    newFunctional->SetGridEnergyWeight( this->m_GridEnergyWeight );
//    newFunctional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    return newFunctional;
    } 
  else
    {
    ImagePairNonrigidRegistrationFunctional *newFunctional = ImagePairNonrigidRegistrationFunctional::Create( this->m_Metric, refVolume, fltVolume, this->m_FloatingImageInterpolation );
    newFunctional->SetAdaptiveFixParameters( this->m_AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->m_AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->m_JacobianConstraintWeight );
    newFunctional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
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
ImagePairNonrigidRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& functional,
  const int idx, const int total ) 
{
  WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->m_Xform );

  float effGridEnergyWeight = this->m_GridEnergyWeight;
  float effJacobianConstraintWeight = this->m_JacobianConstraintWeight;
  float effInverseConsistencyWeight = this->m_InverseConsistencyWeight;

  if ( (this->m_RelaxWeight > 0) && !this->RelaxationStep ) 
    {
    effGridEnergyWeight *= this->m_RelaxWeight;
    effJacobianConstraintWeight *= this->m_RelaxWeight;
    effInverseConsistencyWeight *= this->m_RelaxWeight;
    }

  // handle simple elastic functional
  SmartPointer<ImagePairNonrigidRegistrationFunctional> elasticFunctional = ImagePairNonrigidRegistrationFunctional::SmartPtr::DynamicCastFrom( functional );
  if ( elasticFunctional ) 
    {
    elasticFunctional->SetWarpXform( warpXform );
    elasticFunctional->SetGridEnergyWeight( effGridEnergyWeight );
    elasticFunctional->SetJacobianConstraintWeight( effJacobianConstraintWeight );
    } 
  else 
    {
    // handle inverse-consistent elastic functional
    SmartPointer<ImagePairSymmetricNonrigidRegistrationFunctional> symmetricFunctional = ImagePairSymmetricNonrigidRegistrationFunctional::SmartPtr::DynamicCastFrom( functional );
    if ( symmetricFunctional ) 
      {
      symmetricFunctional->SetWarpXform( warpXform, this->InverseWarpXform );
      symmetricFunctional->SetGridEnergyWeight( effGridEnergyWeight );
      symmetricFunctional->SetJacobianConstraintWeight( effJacobianConstraintWeight );
      symmetricFunctional->SetInverseConsistencyWeight( effInverseConsistencyWeight );
      } 
    else 
      {
      // neither simple nor inverse-consistent functional: something went
      // badly wrong.
      StdErr << "Fatal coding error: Non-elastic functional in ImagePairNonrigidRegistration::EnterResolution.\n";
      abort();
      }
    }
  
  Superclass::EnterResolution( v, functional, idx, total );
}

int 
ImagePairNonrigidRegistration::DoneResolution
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

UniformVolume* 
ImagePairNonrigidRegistration::GetReformattedFloatingImage( Interpolators::InterpolationEnum interpolator )
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
