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

#include "cmtkSplineWarpXform.h"

#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

bool
SplineWarpXform::ApplyInverse
( const Self::SpaceVectorType& v, Self::SpaceVectorType& u, const Types::Coordinate accuracy ) const
{
  u = v;
  return this->ApplyInverseInPlace( u, accuracy );
}

bool
SplineWarpXform::ApplyInverseInPlace 
( Self::SpaceVectorType& v, const Types::Coordinate accuracy ) const
{
  Self::SpaceVectorType u;
  this->FindClosestControlPoint( v, u );
  return this->ApplyInverseInPlaceWithInitial( v, u, accuracy );
}

void
SplineWarpXform::FindClosestControlPoint
( const Self::SpaceVectorType& v, Self::SpaceVectorType& cp ) const
{
  // find closest control point -- we'll go from there.
  Types::Coordinate closestDistance = FLT_MAX;
  Types::Coordinate idx[3];
  for ( int dim = 0; dim < 3; ++dim ) 
    idx[dim] = 0.5 * this->m_Dims[dim];

  for ( Types::Coordinate step = 0.25 * MathUtil::Min( 3, idx ); step > 0.01; step *= 0.5 )
    {
    bool improved = true;
    while ( improved ) 
      {
      improved = false;
      int closestDim = 0, closestDir = 0;
      
      for ( int dim = 0; dim < 3; ++dim ) 
	{
	for ( int dir = -1; dir < 2; dir +=2 )
	  {
	  const Types::Coordinate oldIdx = idx[dim];
	  idx[dim] += dir * step;
	  if ( (idx[dim] > 0) && (idx[dim] <= this->m_Dims[dim]-2) ) 
	    {
	    Self::SpaceVectorType cp;
	    this->GetOriginalControlPointPosition( cp, idx[0], idx[1], idx[2] );
	    this->ApplyInPlace( cp );
	    cp -= v;
	    const Types::Coordinate distance = cp.RootSumOfSquares();
	    if ( distance < closestDistance ) 
	      {
	      closestDistance = distance;
	      closestDim = dim;
	      closestDir = dir;
	      improved = true;
	      } 
	    }
	  idx[dim] = oldIdx;
	  }
	}
      
      if ( improved )
	{
	idx[closestDim] += closestDir * step;
	}
      }
    }
  
  assert( (idx[0] <= this->m_Dims[0]-1) && (idx[1] <= this->m_Dims[1]-1 ) && (idx[2] <= this->m_Dims[2]-1) );
  assert( idx[0] >= 0 && idx[1] >= 0 && idx[2] >= 0 );

  this->GetOriginalControlPointPosition( cp, idx[0], idx[1], idx[2] );
}

bool
SplineWarpXform::ApplyInverseInPlaceWithInitial
( Self::SpaceVectorType& target, const Self::SpaceVectorType& initial, const Types::Coordinate accuracy ) 
  const
{
  Self::SpaceVectorType u( initial );

  // project into domain
  for ( int dim = 0; dim < 3; ++dim )
    {
    u[dim] = std::max<Types::Coordinate>( 0, std::min<Types::Coordinate>( u[dim], this->Domain[dim] ) );
    }

  Self::SpaceVectorType vu( initial ), delta;
  this->ApplyInPlace( vu );
  ((delta = vu) -= target);

  Types::Coordinate error = delta.RootSumOfSquares();

  Types::Coordinate step = 1.0;
  while ( ( error > accuracy) && (step > 0.001) ) 
    {
    // transform difference vector into original coordinate system using inverse Jacobian.
    CoordinateMatrix3x3 J;
    this->GetJacobian( u, J );
    J.Invert3x3();
    J.GetTranspose().Multiply( delta );
    
    // initialize line search
    (vu = u) -= (delta *= step);
    
    // line search along transformed error direction
    if ( !this->InDomain( vu ) ) 
      {
      // project into domain
      for ( int dim = 0; dim < 3; ++dim )
	{
	vu[dim] = std::max<Types::Coordinate>( 0, std::min<Types::Coordinate>( vu[dim], this->Domain[dim] ) );
	}
      }
    
    Self::SpaceVectorType uNext( vu );
    this->ApplyInPlace( vu );
    
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

  target = u;
  return !(error > accuracy);
}

} // namespace cmtk
