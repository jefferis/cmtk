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

#include "cmtkMultiLevelOptimizer.h"

#include <System/cmtkException.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

CallbackResult
MultiLevelOptimizer::Optimize( CoordinateVector& v, const Types::Coordinate, const Types::Coordinate )
{
  if ( ! this->m_Optimizer )
    {
    throw( Exception( "MultiLevelOptimizer.m_Optimizer must be set before calling Optimize().", this ) );
    }

  if ( ! this->m_FunctionalList.size() )
    {
    throw( Exception( "MultiLevelOptimizer must have at least one functional before calling Optimize().", this ) );
    }
  
  // save original value of final functional so we can backtrack later, if
  // need be.
  this->SetFinalValue( (*(this->m_FunctionalList.rbegin()))->m_Functional->EvaluateAt( v ) );
  CoordinateVector vOriginal( v );
  
  // run sequential optimization on all functional with appropriate
  // parameters.
  CallbackResult result = CALLBACK_OK;
  FunctionalListType::iterator fit = this->m_FunctionalList.begin();
  while ( (fit != this->m_FunctionalList.end()) && (result == CALLBACK_OK) )
    {
    this->m_Optimizer->SetFunctional( (*fit)->m_Functional );
    result = this->m_Optimizer->Optimize
      ( v, (*fit)->m_InitialStepSize, (*fit)->m_FinalStepSize );
    ++fit;
    }

  // if we made things worse during the optimization, then
  // backtrack.
  const Self::ReturnType finalFunctionalValue = this->m_Optimizer->GetFinalValue();
  if ( finalFunctionalValue < this->GetFinalValue() )
    {
    v = vOriginal;
    }
  else
    {
    this->SetFinalValue( finalFunctionalValue );
    }
  
  return result;
}

} // namespace cmtk
