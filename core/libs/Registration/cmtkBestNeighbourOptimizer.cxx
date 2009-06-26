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

#include <cmtkBestNeighbourOptimizer.h>

#include <cmtkSearchTrace.h>
#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

CallbackResult 
BestNeighbourOptimizer::Optimize
( CoordinateVector& v, const Self::ParameterType exploration, const Self::ParameterType accuracy )
{
  this->m_LastOptimizeChangedParameters = false;

  int Dim = this->GetSearchSpaceDimension();

  Self::ReturnType optimum = this->Evaluate( v );
  CoordinateVector optimumV(v);

  int optimumDim = -1;
  Types::Coordinate optimumDir = 0;

  const Self::ParameterType real_accuracy = std::min<Self::ParameterType>( exploration, accuracy );
  int numOfSteps = 1+static_cast<int>(log(real_accuracy/exploration)/log(StepFactor));
  Self::ParameterType step = real_accuracy * pow( StepFactor, 1-numOfSteps );
  
  Self::ParameterType *stepScaleVector = Memory::AllocateArray<Self::ParameterType>( Dim );
  for ( int idx=0; idx<Dim; ++idx )
    stepScaleVector[idx] = this->GetParamStep( idx );

  int percentDone = 1, incPercentDone = 99 / numOfSteps;

  SearchTrace<Self::ParameterType> searchTrace ( Dim );

  CallbackResult irq = this->CallbackExecute( v, optimum, percentDone );
  for ( int stepIdx = 0; (stepIdx < numOfSteps) && ( irq == CALLBACK_OK ); ++stepIdx, step *= StepFactor, percentDone += incPercentDone ) 
    {
    char comment[128];
    snprintf( comment, sizeof( comment ), "Setting step size to %4g [mm]", step );
    this->CallbackComment( comment );

    int update = 1;
    while ( update && ( irq == CALLBACK_OK ) ) 
      {
      update = 0;
      
      for ( int dim = 0; dim < Dim; ++dim ) 
	{
	double next;
	
	const Self::ParameterType vOld = v[dim];

	for ( int direction = -1; direction <= 1; direction += 2 )
	  {
	  if ( (irq = this->CallbackExecute( percentDone )) ) break;
	
	  v[dim] = vOld + direction * step * stepScaleVector[dim];
	  if ( !searchTrace.Get( next, dim, step ) )
	    next = this->Evaluate( v );

	  if ( next > optimum ) 
	    {
	    optimum = next;
	    optimumV = v;
	    update = 1;
	    optimumDim = dim;
	    optimumDir = direction * step;
	    }
	  }

	v[dim] = vOld;
	}
      
      if (update) 
	{
	v = optimumV;
	searchTrace.Move( optimumDim, optimumDir );
	irq = this->CallbackExecute( v, optimum, percentDone );
	this->m_LastOptimizeChangedParameters = true;

	// query functional for new parameter steppings if the respective
	// optimizer flag is set.
	if ( this->m_UpdateStepScaleVector )
	  for ( int idx=0; idx<Dim; ++idx )
	    stepScaleVector[idx] = this->GetParamStep( idx );
	}
      }
    }

  this->SetFinalValue( optimum );
  delete[] stepScaleVector;
  return irq;
}

} // namespace cmtk
