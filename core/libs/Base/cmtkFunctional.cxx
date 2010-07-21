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

#include "cmtkFunctional.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

Functional::ReturnType
Functional::EvaluateWithGradient
( CoordinateVector& v, CoordinateVector& g, const Types::Coordinate step )
{ 
  const Self::ReturnType baseValue = this->EvaluateAt( v );
  
  for ( size_t dim = 0; dim < this->VariableParamVectorDim(); ++dim ) 
    {
    const Types::Coordinate stepScale = this->GetParamStep( dim, step );
    if ( stepScale <= 0 ) 
      {
      g[dim] = 0;
      } 
    else
      {
      const Types::Coordinate v0 = v[dim];
      
      v[dim] += stepScale;
      const Self::ReturnType upper = this->EvaluateAt( v );
      
      v[dim] = v0 - stepScale;
      const Self::ReturnType lower = this->EvaluateAt( v );
      
      v[dim] = v0;
      
      if ( (upper > baseValue) || (lower > baseValue) ) 
	{
	g[dim] = upper-lower;
	} 
      else 
	{
	g[dim] = 0;
	}
      }
    }  
  
  return baseValue;
}

} // namespace cmtk
