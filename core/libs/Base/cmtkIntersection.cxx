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

#include "cmtkIntersection.h"

#include <Base/cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

int
Intersection::IntersectX 
( Types::Coordinate& fromFactor, Types::Coordinate& toFactor, 
  const Vector3D& offset, const Vector3D& dX, 
  const Types::Coordinate Size[3], const Types::Coordinate initFromFactor,
  const Types::Coordinate initToFactor, const int lowerClosed, 
  const int upperClosed )
{
  fromFactor = initFromFactor;
  toFactor = initToFactor;
  
  for ( int dim=0; dim<3; ++dim ) 
    {
    if ( dX[dim] > 0 ) 
      {
      fromFactor = std::max( fromFactor, -offset[dim] / dX[dim] );
      toFactor = std::min( toFactor, (Size[dim] - offset[dim]) / dX[dim] );
      } 
    else
      if ( dX[dim] < 0 ) 
	{
	fromFactor = std::max( fromFactor, (Size[dim] - offset[dim]) / dX[dim] );
	toFactor = std::min( toFactor, -offset[dim] / dX[dim] );
	} 
      else
	{
	if ( (offset[dim] < 0) || 
	     ( (offset[dim] == 0) && lowerClosed ) ||
	     (offset[dim] > Size[dim]) ||
	     ( (offset[dim] == Size[dim]) && upperClosed ) ) 
	  {
	  fromFactor = toFactor = 0;
	  return 0;
	  }
	}
    }
  return !( fromFactor > toFactor );
}

int
Intersection::IntersectY 
( Types::Coordinate& fromFactor, Types::Coordinate& toFactor,
  const Vector3D& offset, const Vector3D& dX, 
  const Vector3D& dY, const Types::Coordinate Size[3],
  const Types::Coordinate initFromFactor, const Types::Coordinate initToFactor )
{
  fromFactor = initFromFactor;
  toFactor = initToFactor;

  for ( int dim=0; dim<3; ++dim ) 
    {
    const Types::Coordinate axmin = std::min<Types::Coordinate>( 0, dX[dim] );
    const Types::Coordinate axmax = std::max<Types::Coordinate>( 0, dX[dim] );
    
    if ( dY[dim] > 0 ) 
      {
      fromFactor = std::max( fromFactor, - (offset[dim]+axmax) / dY[dim] );
      toFactor = std::min( toFactor, (Size[dim] - (offset[dim]+axmin)) / dY[dim] );
      } 
    else 
      if ( dY[dim] < 0 ) 
	{
	fromFactor = std::max( fromFactor, (Size[dim] - (offset[dim]+axmin)) / dY[dim] );
	toFactor = std::min( toFactor, - (offset[dim]+axmax) / dY[dim] );
	} 
      else
	{
	if ( (axmax + offset[dim] < 0) || (axmin + offset[dim] > Size[dim]) ) 
	  {
	  fromFactor = toFactor = 0;
	  return 0;
	  }
	}
    }
  return !( fromFactor > toFactor );
}

int
Intersection::IntersectZ 
( Types::Coordinate& fromFactor, Types::Coordinate& toFactor,
  const Vector3D& offset, const Vector3D& dX, 
  const Vector3D& dY, const Vector3D& dZ, 
  const Types::Coordinate Size[3], const Types::Coordinate initFromFactor,
  const Types::Coordinate initToFactor )
{
  fromFactor = initFromFactor;
  toFactor = initToFactor;

  for ( int dim=0; dim<3; ++dim ) 
    {
    const Types::Coordinate axmin = std::min<Types::Coordinate>( 0, dX[dim] );
    const Types::Coordinate axmax = std::max<Types::Coordinate>( 0, dX[dim] );
    const Types::Coordinate bymin = std::min<Types::Coordinate>( 0, dY[dim] );
    const Types::Coordinate bymax = std::max<Types::Coordinate>( 0, dY[dim] );

    if ( dZ[dim] > 0 ) 
      {
      fromFactor = std::max( fromFactor, - (offset[dim]+axmax+bymax) / dZ[dim] );
      toFactor = std::min( toFactor, (Size[dim] - (offset[dim]+axmin+bymin)) / dZ[dim] );
      } 
    else
      if ( dZ[dim] < 0 ) 
	{
	fromFactor = std::max( fromFactor, (Size[dim] - (offset[dim]+axmin+bymin)) / dZ[dim] );
	toFactor = std::min( toFactor, - (offset[dim]+axmax+bymax) / dZ[dim] );
	} 
      else
	{
	if ( (axmax + bymax + offset[dim] < 0) || (axmin + bymin + offset[dim] > Size[dim]) ) 
	  {
	  fromFactor = toFactor = 0;
	  return 0;
	  }
	}
    }
  return !( fromFactor > toFactor );
}

} // namespace cmtk
