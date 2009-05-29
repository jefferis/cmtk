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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkDirectionSetOptimizer.h>

#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

CallbackResult
DirectionSetOptimizer::Optimize
( CoordinateVector& v, const Self::ParameterType exploration, const Self::ParameterType accuracy )
{
  const Self::ParameterType real_accuracy = std::min( exploration, accuracy );

  const Self::ParameterType stepFactor = 0.5;
  Self::ReturnType current = 0;

  const int numOfSteps = 1+static_cast<int>(log(real_accuracy/exploration)/log(stepFactor));
  Self::ParameterType step = real_accuracy * pow( stepFactor, 1-numOfSteps );

  int percentDone = 1, incPercentDone = 99 / numOfSteps;

  CallbackResult irq = CALLBACK_OK;
  for ( int stepIdx = 0; (stepIdx < numOfSteps) && (irq == CALLBACK_OK); ++stepIdx, step *= stepFactor, percentDone += incPercentDone ) 
    {
    char comment[128];
    snprintf( comment, sizeof( comment ), "Setting step size to %4g [mm]", step );
    this->CallbackComment( comment );
    StdErr << comment << "\n";
    
    for ( unsigned int directionIdx = 0; (irq == CALLBACK_OK) && ( directionIdx < this->m_DirectionSet->GetNumberOfDirections()); ++directionIdx ) 
      {
      snprintf( comment, sizeof( comment ), "Selecting direction vector %d", directionIdx );
      this->CallbackComment( comment );
      StdErr << comment << "\n";
      
      Self::ParameterType scaleLength = 1.0 * step;
      
      const CoordinateVector* directionVector = (*this->m_DirectionSet)[directionIdx];
      
      bool update = true;
      while ( update && ( irq == CALLBACK_OK ) ) 
	{
	update = false;
	
	current = this->Evaluate( v );
	irq = this->CallbackExecute( v, current, percentDone );
	
	CoordinateVector vNext( v );
	vNext.Add( *directionVector, scaleLength );
	Self::ReturnType next = this->Evaluate( vNext );
	while (  next > current ) 
	  {
	  if ( ( irq = this->CallbackExecute( percentDone ) ) != CALLBACK_OK )
	    break;
	  current = next;
	  update = true;
	  vNext.Add( *directionVector, scaleLength );
	  next = this->Evaluate( vNext );
	  }
	vNext.Add( *directionVector, -scaleLength );
	if ( update ) v = vNext;
	
	scaleLength *= 0.5;
	
	// Forward-Backward search
	for ( int dirStepIndex = 0; dirStepIndex < numOfSteps; ++dirStepIndex ) 
	  {
	  vNext.Add( *directionVector, scaleLength );
	  Self::ReturnType nextUp = this->Evaluate( vNext );
	  
	  ( vNext = v ).Add( *directionVector, -scaleLength );
	  Self::ReturnType nextDown = this->Evaluate( vNext );
	  
	  if ((nextUp > current) || (nextDown > current)) 
	    {
	    // Here, as we demand strict ">", we prefer smaller steps.
	    if ( nextUp > nextDown ) 
	      {
	      current = nextUp;
	      v.Add( *directionVector, scaleLength );
	      } 
	    else
	      {
	      current = nextDown;
	      v.Add( *directionVector, -scaleLength );
	      }
	    vNext = v;
	    } 
	  
	  scaleLength *= 0.5;
	  
#ifdef DEBUG	      
	  fputs( ".", stderr );
#endif
	  }
	}
#ifdef DEBUG
      fputs( "\n", stderr );
#endif
      
      irq = this->CallbackExecute( v, current, percentDone );
      fprintf( stderr, "%f\n", current );
    }
  }

  return irq;
}

} // namespace cmtk
