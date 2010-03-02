/*
//
//  Copyright 2004-2010 SRI International
//  Copyright 1997-2009 Torsten Rohlfing
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
#include <cmtkTypedArrayFunctionHistogramMatching.h>

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
    m_RepeatMatchFltToRefHistogram( false ),
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
  this->m_FastMode = false;
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

  MatchedLandmarkList::SmartPtr mll( NULL );
  if ( this->m_LandmarkErrorWeight != 0 ) 
    {
    LandmarkList::SmartPtr sourceLandmarks = this->m_ReferenceVolume->m_LandmarkList;
    LandmarkList::SmartPtr targetLandmarks = this->m_FloatingVolume->m_LandmarkList;
    
    if ( sourceLandmarks && targetLandmarks ) 
      {
      this->m_MatchedLandmarks = MatchedLandmarkList::SmartPtr( new MatchedLandmarkList( sourceLandmarks, targetLandmarks ) );
      fprintf( stderr, "Matched %d landmarks.\n", (int)this->m_MatchedLandmarks->size() );
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
  
  if ( this->m_MaxStepSize <= 0 )
    {
    const SplineWarpXform* warp = SplineWarpXform::SmartPtr::DynamicCastFrom( this->m_Xform ); 
    this->m_MaxStepSize = 0.25 * std::max( warp->Spacing[0], std::max( warp->Spacing[1], warp->Spacing[2] ) );
    }

  if ( this->m_CoarsestResolution <= 0 ) 
    this->m_CoarsestResolution = this->m_MaxStepSize;
  
  if ( this->m_UseOriginalData )
    {
    this->m_ParameterStack.push( Self::LevelParameters::SmartPtr( new Self::LevelParameters( -1 ) ) );
    }
  
  for ( ;(currSampling<=this->m_CoarsestResolution); currSampling *= 2 ) 
    {
    this->m_ParameterStack.push( Self::LevelParameters::SmartPtr( new Self::LevelParameters( currSampling ) ) );
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
( const int level, const Superclass::LevelParameters* parameters )
{
  const Self::LevelParameters* levelParameters = dynamic_cast<const Self::LevelParameters*>( parameters );
  if ( ! levelParameters )
    {
    StdErr << "CODING ERROR: wrong RTTI for 'parameters'\n";
    exit( 1 );
    }

  WarpXform::SmartPtr warpXform = WarpXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  if ( ! warpXform )
    {
    StdErr << "CODING ERROR: wrong RTTI for 'this->m_Xform'\n";
    exit( 1 );
    }  

  UniformVolume::SmartPtr referenceVolume( this->m_ReferenceVolume );
  UniformVolume::SmartPtr floatingVolume( this->m_FloatingVolume );
  if ( !level && this->m_MatchFltToRefHistogram )
    {
    floatingVolume = UniformVolume::SmartPtr( floatingVolume->Clone( true /*copyData*/ ) );
    floatingVolume->GetData()->ApplyFunction( TypedArrayFunctionHistogramMatching( *(floatingVolume->GetData()), *(referenceVolume->GetData()) ) );
    }
  else
    {
    if ( this->m_RepeatMatchFltToRefHistogram )
      {
      floatingVolume = UniformVolume::SmartPtr( floatingVolume->Clone( true /*copyData*/ ) );
      UniformVolume::SmartPtr reformat( this->GetReformattedFloatingImage( Interpolators::NEAREST_NEIGHBOR ) );
      floatingVolume->GetData()->ApplyFunction( TypedArrayFunctionHistogramMatching( *(reformat->GetData()), *(referenceVolume->GetData()) ) );
      }
    }
  
  if ( levelParameters->m_Resolution > 0 )
    {
    // for resample if not final, original resolution
    referenceVolume = UniformVolume::SmartPtr( new UniformVolume( *referenceVolume, levelParameters->m_Resolution ) );
    floatingVolume = UniformVolume::SmartPtr( new UniformVolume( *floatingVolume, levelParameters->m_Resolution ) );
    }

  if ( this->m_InverseConsistencyWeight > 0 ) 
    {
    ImagePairSymmetricNonrigidRegistrationFunctional *newFunctional = 
      ImagePairSymmetricNonrigidRegistrationFunctional::Create( this->m_Metric, referenceVolume, floatingVolume, this->m_FloatingImageInterpolation );
    newFunctional->SetInverseConsistencyWeight( this->m_InverseConsistencyWeight );
    newFunctional->SetAdaptiveFixParameters( this->m_AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->m_AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->m_JacobianConstraintWeight );
    newFunctional->SetGridEnergyWeight( this->m_GridEnergyWeight );

    return newFunctional;
    } 
  else
    {
    ImagePairNonrigidRegistrationFunctional *newFunctional = 
      ImagePairNonrigidRegistrationFunctional::Create( this->m_Metric, referenceVolume, floatingVolume, this->m_FloatingImageInterpolation );
    newFunctional->SetAdaptiveFixParameters( this->m_AdaptiveFixParameters );
    newFunctional->SetAdaptiveFixThreshFactor( this->m_AdaptiveFixThreshFactor );
    newFunctional->SetJacobianConstraintWeight( this->m_JacobianConstraintWeight );
    newFunctional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    newFunctional->SetGridEnergyWeight( this->m_GridEnergyWeight );
    if ( this->m_MatchedLandmarks )
      {
      newFunctional->SetLandmarkErrorWeight( this->m_LandmarkErrorWeight );
      newFunctional->SetMatchedLandmarkList( this->m_MatchedLandmarks );
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

  // handle simple nonrigid functional
  SmartPointer<ImagePairNonrigidRegistrationFunctional> nonrigidFunctional = ImagePairNonrigidRegistrationFunctional::SmartPtr::DynamicCastFrom( functional );
  if ( nonrigidFunctional ) 
    {
    nonrigidFunctional->SetWarpXform( warpXform );
    nonrigidFunctional->SetGridEnergyWeight( effGridEnergyWeight );
    nonrigidFunctional->SetJacobianConstraintWeight( effJacobianConstraintWeight );
    } 
  else 
    {
    // handle inverse-consistent nonrigid functional
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
      StdErr << "Fatal coding error: Non-nonrigid functional in ImagePairNonrigidRegistration::EnterResolution.\n";
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
  reformat.SetReferenceVolume( this->m_ReferenceVolume );
  reformat.SetFloatingVolume( this->m_FloatingVolume );

  WarpXform::SmartPtr warpXform( this->GetTransformation() );
  reformat.SetWarpXform( warpXform );

  if ( this->m_ForceOutsideFlag )
    {
    reformat.SetPaddingValue( this->m_ForceOutsideValue );
    }
  
  UniformVolume* result = reformat.PlainReformat();

  if ( this->m_ForceOutsideFlag )
    {
    result->GetData()->ClearPaddingFlag();
    }
  return result;
}


} // namespace cmtk
