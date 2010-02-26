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

#include <cmtkAffineRegistration.h>

#include <cmtkVector.h>

#include <cmtkXform.h>
#include <cmtkAffineXform.h>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>
#include <cmtkFunctional.h>

#include <cmtkVoxelMatchingAffineFunctional.h>
#include <cmtkImagePairAffineRegistrationFunctional.h>
#include <cmtkTypedArrayHistogramMatchingLookup.h>

#include <cmtkOptimizer.h>
#include <cmtkBestNeighbourOptimizer.h>

#include <cmtkReformatVolume.h>

#include <cmtkTimers.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

AffineRegistration::AffineRegistration () :
  m_MatchFltToRefHistogram( false )
{ 
  this->m_InitialAlignCenters = false;
  this->m_NoSwitch = false;
}

AffineRegistration::~AffineRegistration () 
{
}

CallbackResult 
AffineRegistration::InitRegistration ()
{
  CallbackResult result = this->Superclass::InitRegistration();
  if ( result != CALLBACK_OK ) return result;

  if ( this->m_NoSwitch || (this->m_Volume_1->AverageVoxelVolume() >= this->m_Volume_2->AverageVoxelVolume()) ) 
    {
    this->m_ReferenceVolume = this->m_Volume_1;
    this->m_FloatingVolume = this->m_Volume_2;
    SwitchVolumes = false;
    } 
  else
    {
    this->m_ReferenceVolume = this->m_Volume_2;
    this->m_FloatingVolume = this->m_Volume_1;
    SwitchVolumes = true;
    }
  
  if ( this->m_MatchFltToRefHistogram )
    {
    this->GetVolume_2()->GetData()->ApplyLookup( TypedArrayHistogramMatchingLookup( *(this->GetVolume_2()->GetData()), *(this->GetVolume_1()->GetData()) ) );
    }
  
  AffineXform::SmartPtr affineXform;

  if ( this->m_InitialTransformation )
    {
    if ( SwitchVolumes ^ this->m_InitialXformIsInverse )
      {
      affineXform = AffineXform::SmartPtr( this->m_InitialTransformation->MakeInverse() );
      } 
    else
      {
      affineXform = AffineXform::SmartPtr( this->m_InitialTransformation );
      }
  } 
  else 
    {
    affineXform = AffineXform::SmartPtr( new AffineXform );
    }

  if ( this->m_InitialAlignCenters ) 
    {
    Vector3D deltaCenter = ( this->m_FloatingVolume->GetCenterCropRegion() - this->m_ReferenceVolume->GetCenterCropRegion() );
    affineXform->SetXlate( deltaCenter.XYZ );
    }
  
  this->m_Xform = affineXform;
  
  Vector3D center = this->m_ReferenceVolume->GetCenterCropRegion();
  affineXform->ChangeCenter( center.XYZ );

  if ( this->m_UseOriginalData ) 
    {
    Functional::SmartPtr newFunctional( VoxelMatchingAffineFunctional::Create( this->m_Metric, this->m_ReferenceVolume, this->m_FloatingVolume, affineXform ) );
    FunctionalStack.push( newFunctional );
    }
  
  Types::Coordinate currSampling = std::max( this->m_Sampling, 2 * std::min( this->m_ReferenceVolume->GetMinDelta(), this->m_FloatingVolume->GetMinDelta()));
  
  double coarsest = CoarsestResolution;
  if ( coarsest <= 0 ) coarsest = this->m_Exploration;

  UniformVolume::SmartPtr currRef( this->m_ReferenceVolume );
  UniformVolume::SmartPtr currFlt( this->m_FloatingVolume );
  
  for ( ; (currSampling<=coarsest); currSampling *= 2 ) 
    {
    UniformVolume::SmartPtr nextRef( new UniformVolume( *currRef, currSampling ) );
    UniformVolume::SmartPtr nextFlt( new UniformVolume( *currFlt, currSampling ) );
    
    Functional::SmartPtr newFunctional( VoxelMatchingAffineFunctional::Create( this->m_Metric, nextRef, nextFlt, affineXform ) );
    FunctionalStack.push( newFunctional );
    
    currRef = nextRef;
    currFlt = nextFlt;
    }

  this->m_Optimizer = Optimizer::SmartPtr( new BestNeighbourOptimizer( OptimizerStepFactor ) );   
  this->m_Optimizer->SetCallback( this->m_Callback );
  
  // default to rigid transformation
  if ( NumberDOFs.empty() )
    NumberDOFs.push_back( 6 );
  
  // push guard elements
  NumberDOFs.push_back( -1 );
  NumberDOFsFinal.push_back( -1 );
  // intialize iterator.
  NumberDOFsIterator = NumberDOFs.begin();

  return CALLBACK_OK;
}

void
AffineRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, const int level, const int total ) 
{
  if ( *NumberDOFsIterator < 0 )
    {
    if ( (level == total) && (NumberDOFsFinal.size()>1) )
      NumberDOFsIterator = NumberDOFsFinal.begin();
    else
      NumberDOFsIterator = NumberDOFs.begin();
    }

  AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  if ( affineXform ) 
    {
    int numberDOFs = std::min<int>( 12, *NumberDOFsIterator );
    affineXform->SetNumberDOFs( numberDOFs );
    if ( this->m_Callback ) 
      {
      char buffer[64];
      snprintf( buffer, sizeof( buffer ), "Setting Number DOFs to %d.", numberDOFs );
      this->m_Callback->Comment( buffer );
      }
    }
  this->Superclass::EnterResolution( v, f, level, total );
}

int 
AffineRegistration::DoneResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f,
  const int level, const int total )
{
  this->Superclass::DoneResolution( v, f, level, total );
  
  NumberDOFsIterator++;
  return (*NumberDOFsIterator < 0);
}

AffineXform::SmartPtr
AffineRegistration::GetTransformation() const
{
  AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  if ( affineXform && SwitchVolumes ) 
    {
    return affineXform->GetInverse();
    } 
  else 
    {
    return affineXform;
    }
}

UniformVolume* 
AffineRegistration::GetReformattedFloatingImage( Interpolators::InterpolationEnum interpolator )
{
  ReformatVolume reformat;
  reformat.SetInterpolation( interpolator );
  reformat.SetReferenceVolume( this->SwitchVolumes ? this->m_Volume_2 :  this->m_Volume_1 );
  reformat.SetFloatingVolume(  this->SwitchVolumes ? this->m_Volume_1 : this->m_Volume_2 );

  AffineXform::SmartPtr affineXform( this->GetTransformation() );
  reformat.SetAffineXform( affineXform );

  return reformat.PlainReformat();
}

} // namespace cmtk
