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

#include <cmtkBestDirectionOptimizer.h>

#include <cmtkTypes.h>
#include <cmtkConsole.h>
#include <cmtkProgress.h>

#include <algorithm>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

CallbackResult
BestDirectionOptimizer::Optimize
( CoordinateVector& v, const Self::ParameterType exploration, const Self::ParameterType accuracy )
{
  this->m_LastOptimizeChangedParameters = false;
	
  const int Dim = this->GetSearchSpaceDimension();

  const Self::ParameterType real_accuracy = std::min<Self::ParameterType>( exploration, accuracy );
  int numOfSteps = 1+static_cast<int>(log(real_accuracy/exploration)/log(StepFactor));
  Self::ParameterType step = real_accuracy * pow( StepFactor, 1-numOfSteps );

  CoordinateVector directionVector( v.Dim, 0.0 );

  Progress::Begin( 0, numOfSteps, 1, "Multi-resolution optimization" );

  CallbackResult irq = CALLBACK_OK;
  for ( int stepIdx = 0; (stepIdx < numOfSteps) && (irq == CALLBACK_OK); ++stepIdx, step *= StepFactor ) 
    {
    Progress::SetProgress( stepIdx );

    char comment[128];
    snprintf( comment, sizeof( comment ), "Setting step size to %4g [mm]", step );
    this->CallbackComment( comment );
    StdErr.printf( "%s\n", comment );
    
    bool update = true;
    int levelRepeatCounter = this->m_RepeatLevelCount;
    while ( update && ( irq == CALLBACK_OK ) ) 
      {
      update = false;
      
      Self::ReturnType current = this->EvaluateWithGradient( v, directionVector, step );
      irq = this->CallbackExecuteWithData( v, current );

      const Self::ReturnType previous = current;
      
      // Daniel Rueckert is supposedly using Euclid's norm here, but we found this
      // to be less efficient AND accurate. Makes no sense anyway.
      const Self::ParameterType vectorLength = ( this->m_UseMaxNorm ) ? directionVector.MaxNorm() : directionVector.EuclidNorm();
      if ( vectorLength > 0 ) 
	{
	const Self::ParameterType stepLength = step / vectorLength;
	
	// is there a minimum threshold for gradient components? if so,
	// filter out (set to zero) all components below this threshold.
	if ( this->m_DirectionThreshold < 0 ) 
	  {
#pragma omp parallel for
	  for ( int idx=0; idx<Dim; ++idx )
	    directionVector[idx] *= (stepLength * this->GetParamStep(idx) );
	  } 
	else 
	  {
#pragma omp parallel for
	  for ( int idx=0; idx<Dim; ++idx )
	    if ( fabs( directionVector[idx] ) > ( vectorLength * this->m_DirectionThreshold ) ) 
	      {
	      directionVector[idx] *= (stepLength * this->GetParamStep(idx) );
	      } 
	    else 
	      {
	      directionVector[idx] = 0; 
	      }
	  }
	
	CoordinateVector vNext( v );
	vNext += directionVector;
	Self::ReturnType next = this->Evaluate( vNext );
	while (  next > current ) 
	  {
	  if ( ( irq = this->CallbackExecute() ) != CALLBACK_OK )
	    break;
	  current = next;
	  update = true;
	  this->m_LastOptimizeChangedParameters = true;
	  vNext += directionVector;
	  next = this->Evaluate( vNext );
	  }
	vNext -= directionVector;
	if ( update ) v = vNext;
	
	directionVector *= 0.5;
	
	// Forward-Backward search
	for ( int dirStepIndex = 0; dirStepIndex < numOfSteps; ++dirStepIndex ) 
	  {
	  vNext += directionVector;
	  Self::ReturnType nextUp = this->Evaluate( vNext );
	  
	  ( vNext = v ) -= directionVector;
	  Self::ReturnType nextDown = this->Evaluate( vNext );
	  
	  if ((nextUp > current) || (nextDown > current)) 
	    {
	    // Here, as we demand strict ">", we prefer smaller steps.
	    if ( nextUp > nextDown ) 
	      {
	      current = nextUp;
	      v += directionVector;
	      } 
	    else 
	      {
	      current = nextDown;
	      v -= directionVector;
	      }
	    vNext = v;
	    if ( this->m_AggressiveMode )
	      {
	      update = true;
	      this->m_LastOptimizeChangedParameters = true;
	      }
	    } 
	  
	  directionVector *= 0.5;
	  }
	}
      
      irq = this->CallbackExecuteWithData( v, current );
      StdErr.printf( "%f\r", current );

      if ( (fabs(previous-current) / (fabs(previous)+fabs(current)) ) < this->m_DeltaFThreshold )
	update = false;
      
      if ( this->m_AggressiveMode )
	{
	if ( update )
	  {
	  levelRepeatCounter = this->m_RepeatLevelCount;
	  }
	else
	  {
	  --levelRepeatCounter;
	  update = (levelRepeatCounter > 0) && this->m_Functional->Wiggle();
	  }
	}
      }
    }

  Progress::Done();
  
  return irq;
}

} // namespace cmtk
