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

#include <cmtkVoxelRegistration.h>

#include <cmtkVector.h>
#include <cmtkXform.h>
#include <cmtkAffineXform.h>

#include <cmtkFunctional.h>

#include <cmtkTimers.h>
#include <cmtkProgress.h>

#ifdef HAVE_SYS_UTSNAME_H
#  include <sys/utsname.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

VoxelRegistration::VoxelRegistration () 
  : m_PreprocessorRef( "Reference", "ref" ),
    m_PreprocessorFlt( "Floating", "flt" ),
    m_InitialXform( NULL ),
    m_InitialXformIsInverse( false ),
    m_Xform( NULL ),
    m_Optimizer( NULL )
{ 
  this->m_Callback = RegistrationCallback::SmartPtr( new RegistrationCallback() );
  this->m_Protocol = NULL; 

  this->m_Exploration = 8.0;
  this->m_Accuracy = 0.1;
  this->m_Sampling = 1.0;
  CoarsestResolution = -1;
  this->m_UseOriginalData = true;

  this->m_Algorithm = 0;
  UseMaxNorm = true;
  OptimizerStepFactor = 0.5;

  this->m_Metric = 0;
}

VoxelRegistration::~VoxelRegistration () 
{
  if ( this->m_Protocol ) 
    free( this->m_Protocol );
}

CallbackResult
VoxelRegistration::InitRegistration ()
{
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
  
  Types::Coordinate currentExploration = this->m_Exploration;
  CoordinateVector::SmartPtr v( new CoordinateVector() );
  int NumResolutionLevels = FunctionalStack.size();

  Progress::SetTotalSteps( NumResolutionLevels, "Multi-level Registration" );

  int index = 1;
  while ( ! FunctionalStack.empty() && ( irq == CALLBACK_OK ) ) 
    {
    Functional::SmartPtr nextFunctional = FunctionalStack.top();
    FunctionalStack.pop();
    
    // Reference functional as we still need if after the optimization when
    // calling DoneResolution().
    //    nextFunctional->Reference();

    this->m_Optimizer->SetFunctional( nextFunctional );

    int doneResolution = 0;
    while ( ! doneResolution && ( irq == CALLBACK_OK )  ) 
      {
      this->EnterResolution( v, nextFunctional, index, NumResolutionLevels );
      irq = this->ReportProgress( "registration", 1  );
      
      if ( irq == CALLBACK_OK ) 
	{
	Types::Coordinate effectiveAccuracy = (index == NumResolutionLevels) ? std::max<Types::Coordinate>( this->m_Accuracy, currentExploration/1024 ) : this->m_Accuracy;
	
	irq = this->m_Optimizer->Optimize( *v, currentExploration, effectiveAccuracy );
	this->m_Xform->SetParamVector( *v );
	}
      
      doneResolution = this->DoneResolution( v, nextFunctional, index, NumResolutionLevels );
      }
    
    this->m_Optimizer->SetFunctional( Functional::SmartPtr::Null );
    
    currentExploration *= 0.5;

    Progress::SetProgress( index );

    ++index;
    }

  Progress::Done();

  this->OutputResult( v );
  this->ReportProgress( "registration", 100 );
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
    int percentRange = 100 / total;
    this->m_Callback->SetProgressRange( (idx-1) * percentRange, idx * percentRange );
    
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
