/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#include <Registration/cmtkVoxelRegistration.h>

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

VoxelRegistration::VoxelRegistration () 
  : m_Metric( 0 ),
    m_DeltaFThreshold( 0.0 ),
    m_PreprocessorRef( "Reference", "ref" ),
    m_PreprocessorFlt( "Floating", "flt" ),
    m_InitialTransformation( NULL ),
    m_InitialXformIsInverse( false ),
    m_Xform( NULL ),
    m_Optimizer( NULL )
{ 
  this->m_Callback = RegistrationCallback::SmartPtr( new RegistrationCallback() );
  this->m_Protocol = NULL; 

  this->m_Exploration = -1;
  this->m_Accuracy = -1;
  this->m_Sampling = -1;
  this->CoarsestResolution = -1;
  this->m_UseOriginalData = true;

  this->m_Algorithm = 0;
  UseMaxNorm = true;
  OptimizerStepFactor = 0.5;

  this-> SwitchVolumes = false;

  this->TimeStartRegistration = this->TimeStartLevel = this->WalltimeStartRegistration = this->WalltimeStartLevel = this->ThreadTimeStartRegistration = this->ThreadTimeStartLevel = 0.0;
}

VoxelRegistration::~VoxelRegistration () 
{
  free( this->m_Protocol );
}

CallbackResult
VoxelRegistration::InitRegistration ()
{
  if ( this->m_Sampling <= 0 )
    this->m_Sampling = std::max( this->m_Volume_1->GetMaxDelta(), this->m_Volume_2->GetMaxDelta() );
  
  if ( this->m_Exploration <= 0 )
    this->m_Exploration = 8.0 * this->m_Sampling;
  
  if ( this->m_Accuracy <= 0 )
    this->m_Accuracy = this->m_Sampling / 128;
  
  TimeStartLevel = TimeStartRegistration = cmtk::Timers::GetTimeProcess();
  WalltimeStartLevel = WalltimeStartRegistration = cmtk::Timers::GetWalltime();
  ThreadTimeStartLevel = ThreadTimeStartRegistration = cmtk::Timers::GetTimeThread();

  return CALLBACK_OK;
}

CallbackResult
VoxelRegistration::Register ()
{
  CallbackResult irq = this->InitRegistration();
  if ( irq != CALLBACK_OK ) 
    {
    this->DoneRegistration();
    return irq;
    }

  this->m_Optimizer->SetDeltaFThreshold( this->m_DeltaFThreshold );
  
  Types::Coordinate currentExploration = this->m_Exploration;
  CoordinateVector::SmartPtr v( new CoordinateVector() );
  int NumResolutionLevels = FunctionalStack.size();

  Progress::Begin( 0, NumResolutionLevels, 1, "Multi-level Registration" );

  int index = 1;
  while ( ! FunctionalStack.empty() && ( irq == CALLBACK_OK ) ) 
    {
    Functional::SmartPtr nextFunctional = FunctionalStack.top();
    FunctionalStack.pop();
    
    this->m_Optimizer->SetFunctional( nextFunctional );

    int doneResolution = 0;
    while ( ! doneResolution && ( irq == CALLBACK_OK )  ) 
      {
      this->EnterResolution( v, nextFunctional, index, NumResolutionLevels );
      
      if ( irq == CALLBACK_OK ) 
	{
	Types::Coordinate effectiveAccuracy = (index == NumResolutionLevels) ? std::max<Types::Coordinate>( this->m_Accuracy, currentExploration/1024 ) : this->m_Accuracy;
	
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

  this->OutputResult( v, irq );
  this->DoneRegistration( v );
  
  return irq;
}

void
VoxelRegistration::DoneRegistration( const CoordinateVector* v )
{
  if ( v )
    this->m_Xform->SetParamVector( *v );
}

void
VoxelRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, 
  const int idx, const int total ) 
{
  if ( this->m_Callback ) 
    {
    char comment[128];
    snprintf( comment, sizeof( comment ), "Entering resolution level %d out of %d.", idx, total );
    this->m_Callback->Comment( comment );
    }
  
  TimeStartLevel = cmtk::Timers::GetTimeProcess();
  WalltimeStartLevel = cmtk::Timers::GetWalltime();
  ThreadTimeStartLevel = cmtk::Timers::GetTimeThread();
  
  f->GetParamVector( *v );
}

} // namespace cmtk
