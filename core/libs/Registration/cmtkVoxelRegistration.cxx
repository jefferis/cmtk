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

#ifdef HAVE_SYS_UTSNAME_H
#  include <sys/utsname.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

VoxelRegistration::VoxelRegistration () 
  : InitialXform( NULL ),
    InitialXformIsInverse( false ),
    Xform( NULL ),
    Optimizer( NULL )
{ 
  DataClass_1 = DataClass_2 = DATACLASS_GREY;

  Callback = RegistrationCallback::SmartPtr( new RegistrationCallback() );
  Protocol = NULL; 

  Exploration = 8.0;
  Accuracy = 0.1;
  Sampling = 1.0;
  CoarsestResolution = -1;
  UseOriginalData = true;

  Algorithm = 0;
  UseMaxNorm = true;
  OptimizerStepFactor = 0.5;

  Metric = 0;

  // no thresholds (these are FLAGS)
  this->m_ThreshMin1 = this->m_ThreshMin2 = this->m_ThreshMax1 = this->m_ThreshMax2 = 0;
  this->m_ThreshMinValue1 = this->m_ThreshMinValue2 = -FLT_MAX;
  this->m_ThreshMaxValue1 = this->m_ThreshMaxValue2 = FLT_MAX;
}

VoxelRegistration::~VoxelRegistration () 
{
  if ( Protocol ) free( Protocol );
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
  
  Types::Coordinate currentExploration = Exploration;
  CoordinateVector::SmartPtr v( new CoordinateVector() );
  int NumResolutionLevels = FunctionalStack.size();

  int index = 1;
  while ( ! FunctionalStack.empty() && ( irq == CALLBACK_OK ) ) 
    {
    Functional::SmartPtr nextFunctional = FunctionalStack.top();
    FunctionalStack.pop();
    
    // Reference functional as we still need if after the optimization when
    // calling DoneResolution().
    //    nextFunctional->Reference();

    Optimizer->SetFunctional( nextFunctional );

    int doneResolution = 0;
    while ( ! doneResolution && ( irq == CALLBACK_OK )  ) 
      {
      this->EnterResolution( v, nextFunctional, index, NumResolutionLevels );
      irq = this->ReportProgress( "registration", 1  );
      
      if ( irq == CALLBACK_OK ) 
	{
	Types::Coordinate effectiveAccuracy = (index == NumResolutionLevels) ? std::max<Types::Coordinate>( Accuracy, currentExploration/1024 ) : Accuracy;
	
	irq = Optimizer->Optimize( *v, currentExploration, effectiveAccuracy );
	Xform->SetParamVector( *v );
	}
      
      doneResolution = this->DoneResolution( v, nextFunctional, index, NumResolutionLevels );
      }
    
    Optimizer->SetFunctional( Functional::SmartPtr::Null );
    
    currentExploration *= 0.5;
    ++index;
    }
  this->OutputResult( v );
  this->ReportProgress( "registration", 100 );
  this->DoneRegistration( v );
  
  return irq;
}

void
VoxelRegistration::DoneRegistration( const CoordinateVector* v )
{
  if ( v )
    Xform->SetParamVector( *v );
}

void
VoxelRegistration::EnterResolution
( CoordinateVector::SmartPtr& v, Functional::SmartPtr& f, 
  const int idx, const int total ) 
{
  if ( Callback ) 
    {
    int percentRange = 100 / total;
    Callback->SetProgressRange( (idx-1) * percentRange, idx * percentRange );
    
    char comment[128];
    snprintf( comment, sizeof( comment ), "Entering resolution level %d out of %d.", idx, total );
    Callback->Comment( comment );
    }
  
  TimeStartLevel = cmtk::Timers::GetTimeProcess();
  WalltimeStartLevel = cmtk::Timers::GetWalltime();
  ThreadTimeStartLevel = cmtk::Timers::GetTimeThread();
  
  f->GetParamVector( *v );
}

} // namespace cmtk
