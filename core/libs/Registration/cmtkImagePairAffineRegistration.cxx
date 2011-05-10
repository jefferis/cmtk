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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include "cmtkImagePairAffineRegistration.h"

#include <System/cmtkTimers.h>

#include <Base/cmtkVector.h>
#include <Base/cmtkXform.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkTypedArrayFunctionHistogramMatching.h>
#include <Base/cmtkVolume.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkFunctional.h>

#include <Registration/cmtkImagePairAffineRegistrationFunctional.h>
#include <Registration/cmtkImagePairSymmetricAffineRegistrationFunctional.h>
#include <Registration/cmtkOptimizer.h>
#include <Registration/cmtkBestNeighbourOptimizer.h>
#include <Registration/cmtkReformatVolume.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairAffineRegistration::ImagePairAffineRegistration () :
  m_Initializer( MakeInitialAffineTransformation::FOV ),
  m_SymmetricFwdBwd( false ),
  m_MatchFltToRefHistogram( false )
{ 
}

ImagePairAffineRegistration::~ImagePairAffineRegistration () 
{
}

CallbackResult 
ImagePairAffineRegistration::InitRegistration ()
{
  CallbackResult result = this->Superclass::InitRegistration();
  if ( result != CALLBACK_OK ) return result;

  this->m_ReferenceVolume = this->m_Volume_1;
  this->m_FloatingVolume = this->m_Volume_2;
  
  if ( this->m_MatchFltToRefHistogram )
    {
    this->GetVolume_2()->GetData()->ApplyFunctionObject( TypedArrayFunctionHistogramMatching( *(this->GetVolume_2()->GetData()), *(this->GetVolume_1()->GetData()) ) );
    }
  
  AffineXform::SmartPtr affineXform;

  if ( this->m_InitialTransformation )
    {
    if ( this->m_InitialXformIsInverse )
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
    affineXform = AffineXform::SmartPtr( MakeInitialAffineTransformation::Create( *this->m_ReferenceVolume, *this->m_FloatingVolume, this->m_Initializer ) );
    }
  
  this->m_Xform = affineXform;
  
  Vector3D center = this->m_ReferenceVolume->GetCenterCropRegion();
  affineXform->ChangeCenter( center );

  if ( this->m_UseOriginalData ) 
    {
    this->m_ParameterStack.push( Self::LevelParameters::SmartPtr( new Self::LevelParameters( -1 ) ) );
    }
  
  Types::Coordinate currSampling = std::max( this->m_Sampling, 2 * std::min( this->m_ReferenceVolume->GetMinDelta(), this->m_FloatingVolume->GetMinDelta()));
  
  double coarsest = this->m_CoarsestResolution;
  if ( coarsest <= 0 ) coarsest = this->m_MaxStepSize;

  for ( ; (currSampling<=coarsest); currSampling *= 2 ) 
    {
    this->m_ParameterStack.push( Self::LevelParameters::SmartPtr( new Self::LevelParameters( currSampling ) ) );
    }
  
  this->m_Optimizer = Optimizer::SmartPtr( new BestNeighbourOptimizer( this->m_OptimizerStepFactor ) );   
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

Functional* 
ImagePairAffineRegistration
::MakeFunctional( const int /*level*/, const Superclass::LevelParameters* parameters )
{
  const Self::LevelParameters* levelParameters = dynamic_cast<const Self::LevelParameters*>( parameters );
  if ( ! levelParameters )
    {
    StdErr << "CODING ERROR: wrong RTTI for 'parameters'\n";
    exit( 1 );
    }

  AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  if ( ! affineXform )
    {
    StdErr << "CODING ERROR: wrong RTTI for 'this->m_Xform'\n";
    exit( 1 );
    }
  
  UniformVolume::SmartPtr nextRef, nextFlt;
  if ( levelParameters->m_Resolution > 0 )
    {
    nextRef = UniformVolume::SmartPtr( new UniformVolume( *this->m_ReferenceVolume, levelParameters->m_Resolution ) );
    nextFlt = UniformVolume::SmartPtr( new UniformVolume( *this->m_FloatingVolume, levelParameters->m_Resolution ) );
    }
  else
    {
    // for final, original resolution just take input volumes.
    nextRef = this->m_ReferenceVolume;
    nextFlt = this->m_FloatingVolume;
    }

  if ( this->m_SymmetricFwdBwd )
    {
    ImagePairSymmetricAffineRegistrationFunctional *functional = ImagePairSymmetricAffineRegistrationFunctional::Create( this->m_Metric, nextRef, nextFlt, this->m_FloatingImageInterpolation, affineXform );
    functional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    return functional;
    }
  else
    {
    ImagePairAffineRegistrationFunctional *functional = ImagePairAffineRegistrationFunctional::Create( this->m_Metric, nextRef, nextFlt, this->m_FloatingImageInterpolation, affineXform );
    functional->SetForceOutside( this->m_ForceOutsideFlag, this->m_ForceOutsideValue );
    return functional;
    }
}

void
ImagePairAffineRegistration::EnterResolution
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
    affineXform->SetNumberDOFs( *NumberDOFsIterator );
    if ( this->m_Callback ) 
      {
      char buffer[64];
      snprintf( buffer, sizeof( buffer ), "Setting Number DOFs to %d.",  *NumberDOFsIterator );
      this->m_Callback->Comment( buffer );
      }
    }
  this->Superclass::EnterResolution( v, f, level, total );
}

int 
ImagePairAffineRegistration::DoneResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f,
  const int level, const int total )
{
  this->Superclass::DoneResolution( v, f, level, total );
  
  NumberDOFsIterator++;
  return (*NumberDOFsIterator < 0);
}

AffineXform::SmartPtr
ImagePairAffineRegistration::GetTransformation() const
{
  AffineXform::SmartPtr affineXform = AffineXform::SmartPtr::DynamicCastFrom( this->m_Xform );
  return affineXform;
}

const UniformVolume::SmartPtr
ImagePairAffineRegistration::GetReformattedFloatingImage( Interpolators::InterpolationEnum interpolator ) const
{
  ReformatVolume reformat;
  reformat.SetInterpolation( interpolator );
  reformat.SetReferenceVolume( this->m_Volume_1 );
  reformat.SetFloatingVolume( this->m_Volume_2 );

  AffineXform::SmartPtr affineXform( this->GetTransformation() );
  reformat.SetAffineXform( affineXform );

  return reformat.PlainReformat();
}

} // namespace cmtk
