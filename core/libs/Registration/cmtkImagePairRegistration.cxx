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

#include "cmtkImagePairRegistration.h"

#include <Base/cmtkVector.h>
#include <Base/cmtkXform.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkFunctional.h>

#include <System/cmtkTimers.h>
#include <System/cmtkProgress.h>

#ifdef HAVE_SYS_UTSNAME_H
#  include <sys/utsname.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

ImagePairRegistration::ImagePairRegistration () 
  : m_Metric( 0 ),
    m_FloatingImageInterpolation( Interpolators::DEFAULT ),
    m_AutoMultiLevels( 0 ),
    m_MaxStepSize( -1 ),
    m_MinStepSize( -1 ),
    m_DeltaFThreshold( 0.0 ),
    m_Sampling( -1 ),
    m_ForceOutsideFlag( false ),
    m_ForceOutsideValue( 0.0 ),
    m_PreprocessorRef( "Reference", "ref" ),
    m_PreprocessorFlt( "Floating", "flt" ),
    m_InitialTransformation( NULL ),
    m_InitialXformIsInverse( false ),
    m_Xform( NULL ),
    m_Optimizer( NULL )
{ 
  this->m_Callback = RegistrationCallback::SmartPtr( new RegistrationCallback() );

  this->m_Sampling = -1;
  this->m_CoarsestResolution = -1;
  this->m_UseOriginalData = true;

  this->m_Algorithm = 0;
  this->m_UseMaxNorm = true;
  this->m_OptimizerStepFactor = 0.5;
}

CallbackResult
ImagePairRegistration::InitRegistration ()
{
  if ( this->m_AutoMultiLevels > 0 )
    {
    const Types::Coordinate minDelta = std::min( this->m_Volume_1->GetMinDelta(), this->m_Volume_2->GetMinDelta() );
    const Types::Coordinate maxDelta = std::max( this->m_Volume_1->GetMaxDelta(), this->m_Volume_2->GetMaxDelta() );

    this->m_MinStepSize = 0.1 * minDelta;
    this->m_Sampling = maxDelta;
    this->m_MaxStepSize = maxDelta * (1<<(this->m_AutoMultiLevels-1));
    }
  
  if ( this->m_Sampling <= 0 )
    this->m_Sampling = std::max( this->m_Volume_1->GetMaxDelta(), this->m_Volume_2->GetMaxDelta() );
  
  if ( this->m_MaxStepSize <= 0 )
    this->m_MaxStepSize = 8.0 * this->m_Sampling;
  
  if ( this->m_MinStepSize <= 0 )
    this->m_MinStepSize = this->m_Sampling / 128;
  
  this->m_TimeStartLevel = this->m_TimeStartRegistration = cmtk::Timers::GetTimeProcess();
  this->m_WalltimeStartLevel = this->m_WalltimeStartRegistration = cmtk::Timers::GetWalltime();
  this->m_ThreadTimeStartLevel = this->m_ThreadTimeStartRegistration = cmtk::Timers::GetTimeThread();

  return CALLBACK_OK;
}

CallbackResult
ImagePairRegistration::Register ()
{
  CallbackResult irq = this->InitRegistration();
  if ( irq != CALLBACK_OK ) 
    {
    this->DoneRegistration();
    return irq;
    }

  this->m_Optimizer->SetDeltaFThreshold( this->m_DeltaFThreshold );
  
  Types::Coordinate currentExploration = this->m_MaxStepSize;
  CoordinateVector::SmartPtr v( new CoordinateVector() );
  const size_t NumResolutionLevels = this->m_ParameterStack.size();
  
  Progress::Begin( 0, NumResolutionLevels, 1, "Multi-level Registration" );

  unsigned int index = 1;
  while ( ! this->m_ParameterStack.empty() && ( irq == CALLBACK_OK ) ) 
    {
    Functional::SmartPtr nextFunctional( this->MakeFunctional( index-1, this->m_ParameterStack.top() ) );
    this->m_ParameterStack.pop();
    
    // Reference functional as we still need if after the optimization when
    // calling DoneResolution().
    //    nextFunctional->Reference();
    
    this->m_Optimizer->SetFunctional( nextFunctional );

    int doneResolution = 0;
    while ( ! doneResolution && ( irq == CALLBACK_OK )  ) 
      {
      this->EnterResolution( v, nextFunctional, index, NumResolutionLevels );
      
      if ( irq == CALLBACK_OK ) 
	{
	Types::Coordinate effectiveAccuracy = (index == NumResolutionLevels) ? std::max<Types::Coordinate>( this->m_MinStepSize, currentExploration/1024 ) : this->m_MinStepSize;
	
	irq = this->m_Optimizer->Optimize( *v, currentExploration, effectiveAccuracy );
	this->m_Xform->SetParamVector( *v );
	}
      
      doneResolution = this->DoneResolution( v, nextFunctional, index, NumResolutionLevels );
      }
    
    this->m_Optimizer->SetFunctional( Functional::SmartPtr::Null() );
    
    currentExploration *= 0.5;

    Progress::SetProgress( index );

    ++index;
    }

  Progress::Done();

  this->OutputResult( v );
  this->DoneRegistration( v );
  
  return irq;
}

void
ImagePairRegistration::DoneRegistration( const CoordinateVector* v )
{
  if ( v )
    this->m_Xform->SetParamVector( *v );
}

void
ImagePairRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, 
  const int idx, const int total ) 
{
  if ( this->m_Callback ) 
    {
    char comment[128];
    snprintf( comment, sizeof( comment ), "Entering resolution level %d out of %d.", idx, total );
    this->m_Callback->Comment( comment );
    }
  
  this->m_TimeStartLevel = cmtk::Timers::GetTimeProcess();
  this->m_WalltimeStartLevel = cmtk::Timers::GetWalltime();
  this->m_ThreadTimeStartLevel = cmtk::Timers::GetTimeThread();
  
  f->GetParamVector( *v );
}

} // namespace cmtk
