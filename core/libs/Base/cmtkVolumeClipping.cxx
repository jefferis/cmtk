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

#include "cmtkVolumeClipping.h"

#include <Base/cmtkMathUtil.h>

#include <string.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void 
VolumeClipping::SetClippingBoundaries
(  const UniformVolume::CoordinateRegionType& region )
{
  this->m_ClippingRegion = region;
}

int VolumeClipping::ClipX 
( Types::Coordinate& fromFactor, Types::Coordinate& toFactor,
  const Vector3D& offset, 
  const Types::Coordinate initFromFactor, const Types::Coordinate initToFactor,
  const bool lowerClosed, const bool upperClosed ) const
{
  fromFactor = initFromFactor;
  toFactor = initToFactor;
  
  for ( int dim=0; dim<3; ++dim ) 
    {
    if ( DeltaX[dim] > 0 ) 
      {
      fromFactor = std::max( fromFactor, (this->m_ClippingRegion.From()[dim] - offset[dim]) / DeltaX[dim] );
      toFactor = std::min( toFactor, (this->m_ClippingRegion.To()[dim] - offset[dim]) / DeltaX[dim] );
      } 
    else
      if ( DeltaX[dim] < 0 ) 
	{
	fromFactor = std::max( fromFactor, (this->m_ClippingRegion.To()[dim] - offset[dim]) / DeltaX[dim] );
	toFactor = std::min( toFactor, (this->m_ClippingRegion.From()[dim]-offset[dim]) / DeltaX[dim] );
	} 
      else
	{
	if ( (offset[dim] < this->m_ClippingRegion.From()[dim]) || 
	     ( (offset[dim] == this->m_ClippingRegion.From()[dim] ) && lowerClosed ) ||
	     (offset[dim] > this->m_ClippingRegion.To()[dim]) ||
	     ( (offset[dim] == this->m_ClippingRegion.To()[dim]) && upperClosed ) ) {
	fromFactor = toFactor = 0;
	return 0;
	}
	}
    }
  return !( fromFactor > toFactor );
}

int VolumeClipping::ClipY
( Types::Coordinate& fromFactor, Types::Coordinate& toFactor,
  const Vector3D& offset, 
  const Types::Coordinate initFromFactor, const Types::Coordinate initToFactor ) const
{
  fromFactor = initFromFactor;
  toFactor = initToFactor;

  for ( int dim=0; dim<3; ++dim ) 
    {
#ifdef _MSC_VER
    const Types::Coordinate axmin = offset[dim] + std::min( (Types::Coordinate) 0, DeltaX[dim] );
    const Types::Coordinate axmax = offset[dim] + std::max( (Types::Coordinate) 0, DeltaX[dim] );
#else
    const Types::Coordinate axmin = offset[dim] + std::min<Types::Coordinate>( 0, DeltaX[dim] );
    const Types::Coordinate axmax = offset[dim]+ std::max<Types::Coordinate>( 0, DeltaX[dim] );
#endif

    if ( DeltaY[dim] > 0 ) 
      {
      fromFactor = std::max( fromFactor, ( this->m_ClippingRegion.From()[dim] - axmax) / DeltaY[dim] );
      toFactor = std::min( toFactor, ( this->m_ClippingRegion.To()[dim] - axmin) / DeltaY[dim] );
      } 
    else
      if ( DeltaY[dim] < 0 ) 
	{
	fromFactor = std::max( fromFactor, (this->m_ClippingRegion.To()[dim] - axmin) / DeltaY[dim] );
	toFactor = std::min( toFactor, (this->m_ClippingRegion.From()[dim] - axmax) / DeltaY[dim] );
	} 
      else
	{
	if ( (axmax < this->m_ClippingRegion.From()[dim]) || (axmin > this->m_ClippingRegion.To()[dim]) ) 
	  {
	  fromFactor = toFactor = 0;
	  return 0;
	  }
	}
    }
  return !( fromFactor > toFactor );
}

int VolumeClipping::ClipZ
( Types::Coordinate& fromFactor, Types::Coordinate& toFactor, 
  const Vector3D& offset, 
  const Types::Coordinate initFromFactor, const Types::Coordinate initToFactor ) const
{
  fromFactor = initFromFactor;
  toFactor = initToFactor;

  for ( int dim=0; dim<3; ++dim ) 
    {
#ifdef _MSC_VER
    const Types::Coordinate axbymin = offset[dim] + std::min( (Types::Coordinate) 0, DeltaX[dim] ) + std::min( (Types::Coordinate) 0, DeltaY[dim] );
    const Types::Coordinate axbymax = offset[dim] + std::max( (Types::Coordinate) 0, DeltaX[dim] ) + std::max( (Types::Coordinate) 0, DeltaY[dim] );
#else
    const Types::Coordinate axbymin = offset[dim] + std::min<Types::Coordinate>( 0, DeltaX[dim] ) + std::min<Types::Coordinate>( 0, DeltaY[dim] );
    const Types::Coordinate axbymax = offset[dim] + std::max<Types::Coordinate>( 0, DeltaX[dim] ) + std::max<Types::Coordinate>( 0, DeltaY[dim] );
#endif

    if ( DeltaZ[dim] > 0 ) 
      {
      fromFactor = std::max( fromFactor, (this->m_ClippingRegion.From()[dim] - axbymax) / DeltaZ[dim] );
      toFactor = std::min( toFactor, (this->m_ClippingRegion.To()[dim] - axbymin) / DeltaZ[dim] );
      } 
    else
      if ( DeltaZ[dim] < 0 ) 
	{
	fromFactor = std::max( fromFactor, (this->m_ClippingRegion.To()[dim] - axbymin) / DeltaZ[dim] );
	toFactor = std::min( toFactor, (this->m_ClippingRegion.From()[dim] - axbymax) / DeltaZ[dim] );
	} 
      else
	{
	if ( (axbymax < this->m_ClippingRegion.From()[dim]) || (axbymin > this->m_ClippingRegion.To()[dim]) ) 
	  {
	  fromFactor = toFactor = 0;
	  return 0;
	  }
	}
    }
  return !( fromFactor > toFactor );
}

}  // namespace cmtk
