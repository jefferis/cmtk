/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkUniformVolumePainter.h"

void
cmtk::UniformVolumePainter::DrawSphere
( const UniformVolume::CoordinateVectorType& center, const Types::Coordinate radius, const Types::DataItem value )
{
  UniformVolume& volume = *(this->m_Volume);

  UniformVolume::CoordinateVectorType centerAbsolute( center );
  Types::Coordinate radiusAbsolute[3] = { radius, radius, radius };
  
  switch ( this->m_CoordinateMode )
    {
    default:
    case Self::COORDINATES_ABSOLUTE:
      // nothing to do - already absolute
      break;
    case Self::COORDINATES_RELATIVE:
      for ( int dim = 0; dim < 3; ++dim )
	{
	(centerAbsolute[dim] *= volume.Size[dim]) += volume.m_Offset[dim];
	radiusAbsolute[dim] *= volume.Size[dim];
	}
      break;
    case Self::COORDINATES_INDEXED:
      for ( int dim = 0; dim < 3; ++dim )
	{
	(centerAbsolute[dim] *= volume.m_Delta[dim]) += volume.m_Offset[dim];
	radiusAbsolute[dim] *= volume.m_Delta[dim];
	}
      break;
    }
  
  DataGrid::RegionType region = volume.GetWholeImageRegion();
  for ( int dim = 0; dim < 3; ++dim )
    {
    region.From()[dim] = std::max( region.From()[dim], static_cast<int>( (centerAbsolute[dim] - radiusAbsolute[dim] - volume.m_Offset[dim]) / volume.m_Delta[dim] ) - 1 );
    region.To()[dim] = std::min( region.To()[dim], static_cast<int>( (centerAbsolute[dim] + radiusAbsolute[dim] - volume.m_Offset[dim]) / volume.m_Delta[dim] ) + 1 );
    }

  for ( int k = region.From()[2]; k < region.To()[2]; ++k )
    {
    const Types::Coordinate Z = volume.GetPlaneCoord( 2, k );
    for ( int j = region.From()[1]; j < region.To()[1]; ++j )
      {
      const Types::Coordinate Y = volume.GetPlaneCoord( 1, j );
      
      size_t offset = region.From()[0] + volume.m_Dims[0] * ( j + volume.m_Dims[1] * k );
      for ( int i = region.From()[0]; i < region.To()[0]; ++i, ++offset )
	{
	const Types::Coordinate X = volume.GetPlaneCoord( 0, i );
	
	UniformVolume::CoordinateVectorType v = FixedVectorStaticInitializer<3,Types::Coordinate>::Init( X, Y, Z );
	v -= centerAbsolute;

	for ( int dim = 0; dim < 3; ++dim )
	  {
	  v[dim] /= radiusAbsolute[dim];
	  }

	if ( v.RootSumOfSquares() <= 1.0 )
	  volume.SetDataAt( value, offset );
	}
      }
    }
}

void
cmtk::UniformVolumePainter::DrawBox
( const UniformVolume::CoordinateVectorType& boxFrom, const UniformVolume::CoordinateVectorType& boxTo, const Types::DataItem value )
{
  UniformVolume& volume = *(this->m_Volume);

  int indexFrom[3], indexTo[3];

  switch ( this->m_CoordinateMode )
    {
    default:
    case Self::COORDINATES_ABSOLUTE:
      for ( int dim = 0; dim < 3; ++dim )
	{
	indexFrom[dim] = static_cast<int>( MathUtil::Round( boxFrom[dim] / volume.m_Delta[dim] ) );
	indexTo[dim] = static_cast<int>( MathUtil::Round( boxTo[dim] / volume.m_Delta[dim] ) );
	}
      break;
    case Self::COORDINATES_RELATIVE:
      for ( int dim = 0; dim < 3; ++dim )
	{
	indexFrom[dim] = static_cast<int>( MathUtil::Round( boxFrom[dim] * volume.Size[dim] / volume.m_Delta[dim] ) );
	indexTo[dim] = static_cast<int>( MathUtil::Round( boxTo[dim] * volume.Size[dim] / volume.m_Delta[dim] ) );
	}
      break;
    case Self::COORDINATES_INDEXED:
      // nothing to do - already indexed
      for ( int dim = 0; dim < 3; ++dim )
	{
	indexFrom[dim] = static_cast<int>( boxFrom[dim] + 0.5 );
	indexTo[dim] = static_cast<int>( boxTo[dim] + 0.5 );
	}
      break;
    }

  // make sure boundaries are in correct order and in valid range
  for ( int dim = 0; dim < 3; ++dim )
    {
    if ( indexFrom[dim] > indexTo[dim] )
      std::swap( indexFrom[dim], indexTo[dim] );

    indexFrom[dim] = std::max( 0, std::min( volume.m_Dims[dim]-1, indexFrom[dim] ) );
    indexTo[dim] = std::max( 0, std::min( volume.m_Dims[dim]-1, indexTo[dim] ) );
    }
  
  for ( int k = indexFrom[2]; k <= indexTo[2]; ++k )
    {
    for ( int j = indexFrom[1]; j <= indexTo[1]; ++j )
      {
      for ( int i = indexFrom[0]; i <= indexTo[0]; ++i )
	{
	volume.SetDataAt( value, volume.GetOffsetFromIndex( i, j, k ) );
	}
      }
    }
}
