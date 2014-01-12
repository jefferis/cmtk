/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
//  $Revision: 4929 $
//
//  $LastChangedDate: 2013-10-04 10:18:35 -0700 (Fri, 04 Oct 2013) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include "cmtkXform.h"

bool
cmtk::Xform::ApplyInverseWithInitial
( const Self::SpaceVectorType& target, Self::SpaceVectorType& source, const Self::SpaceVectorType& initial, const Types::Coordinate accuracy ) 
  const
{
  Self::SpaceVectorType u( initial );
  this->ProjectToDomain( u );

  Self::SpaceVectorType vu( this->Apply( initial ) ), delta;
  ((delta = vu) -= target);

  Types::Coordinate error = delta.RootSumOfSquares();

  Types::Coordinate step = 1.0;
  while ( ( error > accuracy) && (step > 0.001) ) 
    {
    // transform difference vector into original coordinate system using inverse Jacobian.
    delta *= this->GetJacobian( u ).GetInverse().GetTranspose();
    
    // initialize line search
    (vu = u) -= (delta *= step);

    // project back into transformation domain, if necessary
    this->ProjectToDomain( vu );
    
    // line search along transformed error direction
    Self::SpaceVectorType uNext( vu );
    vu = this->Apply( vu );
    
    (delta = vu) -= target;
    if ( error > delta.RootSumOfSquares() ) 
      {
      error = delta.RootSumOfSquares();
      u = uNext;
      } 
    else
      {
      step *= 0.5;
      }
    }

  source = u;
  return !(error > accuracy);
}
