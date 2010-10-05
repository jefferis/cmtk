/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkIterativeDirectionOptimizer.h"

#include <System/cmtkProgress.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

CallbackResult 
IterativeDirectionOptimizer::Optimize
( CoordinateVector& v, const Types::Coordinate exploration, const Types::Coordinate accuracy )
{
  int Dim = this->GetSearchSpaceDimension();

  Types::Coordinate optimum = this->Evaluate( v );
  CoordinateVector optimumV(v);

  const double real_accuracy = std::min<Types::Coordinate>( exploration, accuracy );
  int numOfSteps = 1+static_cast<int>(log(real_accuracy/exploration)/log(StepFactor));
  Types::Coordinate step = real_accuracy * pow( StepFactor, 1-numOfSteps );

  std::vector<Types::Coordinate> stepScaleVector( Dim );
  for ( int idx=0; idx<Dim; ++idx )
    stepScaleVector[idx] = this->GetParamStep( idx );

  Progress::Begin( 0, numOfSteps, 1, "Multi-resolution optimization" );

  CallbackResult irq = this->CallbackExecuteWithData( v, optimum );
  for ( int stepIdx = 0; (stepIdx < numOfSteps) && ( irq == CALLBACK_OK ); ++stepIdx, step *= StepFactor ) 
    {
    Progress::SetProgress( stepIdx );

    char comment[128];
    snprintf( comment, sizeof( comment ), "Setting step size to %4g [mm]", step );
    this->CallbackComment( comment );

    bool update = true;
    while ( update && ( irq == CALLBACK_OK ) ) 
      {
      update = false;
      
      for ( int dim = 0; (dim < Dim) && ( irq == CALLBACK_OK ); ++dim ) 
	{
	bool updateThisDim = true;
	while ( updateThisDim && ( irq == CALLBACK_OK ) ) 
	  {
	  updateThisDim = false;
	  
	  Types::Coordinate vOld = v[dim];
	  
	  if ( (irq = this->CallbackExecute()) ) break;
	  
	  v[dim] += step * stepScaleVector[dim];
      const Self::ReturnType fUpper = this->Evaluate( v );
	  
	  Types::Coordinate optimumStep = 0;
	  
	  if ( fUpper > optimum ) 
	    {
	    optimum = fUpper;
	    optimumV = v;
	    updateThisDim = true;
	    optimumStep = step;
	    }
	  
	  if ( (irq = this->CallbackExecute()) ) break;
	  
	  v[dim] = vOld - (step * stepScaleVector[dim]);
      const Self::ReturnType fLower = this->Evaluate( v );
	  
	  if ( fLower > optimum ) 
	    {
	    optimum = fLower;
	    optimumV = v;
	    updateThisDim = true;
	    optimumStep = -step;
	    }
	  
	  if ( updateThisDim ) 
	    {
	    update = true;
	    while ( updateThisDim ) 
	      {
	      updateThisDim = false;
	      vOld = v[dim];
	      v[dim] += optimumStep * stepScaleVector[dim];
          const Self::ReturnType f = this->Evaluate( v );
	      
	      if ( f > optimum ) 
		{
		optimum = f;
		optimumV = v;
		updateThisDim = true;
		}
	      
	      if ( (irq = this->CallbackExecute()) ) break;
	      }
	    updateThisDim = true;
	    }
	  
	  v[dim] = vOld;
	  }
	
	if ( update ) 
	  {
	  v = optimumV;
	  irq = this->CallbackExecuteWithData( v, optimum );
	  // query functional for new parameter steppings if the respective
	  // optimzier flag is set.
	  if ( this->m_UpdateStepScaleVector )
	    for ( int idx=0; idx<Dim; ++idx )
	      stepScaleVector[idx] = this->GetParamStep( idx );
	  }
	}
      }
    }

  Progress::Done();
  
  return irq;
}

} // namespace cmtk
