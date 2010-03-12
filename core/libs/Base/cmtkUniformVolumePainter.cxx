/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkUniformVolumePainter.h>

void
cmtk::UniformVolumePainter::DrawSphere
( const Vector3D& center, const Types::Coordinate radius, const Types::DataItem value )
{
  Vector3D centerAbsolute( center );
  Types::Coordinate radiusAbsolute[3] = { radius, radius, radius };
  
  switch ( this->m_CoordinateMode )
    {
    default:
    case Self::ABSOLUTE:
      // nothing to do - already absolute
      break;
    case Self::RELATIVE:
      for ( int dim = 0; dim < 3; ++dim )
	{
	centerAbsolute[dim] *= this->m_Volume->Size[dim];
	radiusAbsolute[dim] *= this->m_Volume->Size[dim];
	}
      break;
    case Self::INDEXED:
      for ( int dim = 0; dim < 3; ++dim )
	{
	centerAbsolute[dim] *= this->m_Volume->m_Delta[dim];
	radiusAbsolute[dim] *= this->m_Volume->m_Delta[dim];
	}
      break;
    }
  
  size_t offset = 0;
  for ( int k = 0; k < this->m_Volume->m_Dims[2]; ++k )
    {
    const Types::Coordinate Z = this->m_Volume->GetPlaneCoord( 2, k );
    for ( int j = 0; j < this->m_Volume->m_Dims[1]; ++j )
      {
    const Types::Coordinate Y = this->m_Volume->GetPlaneCoord( 1, j );
      for ( int i = 0; i < this->m_Volume->m_Dims[0]; ++i, ++offset )
	{
	const Types::Coordinate X = this->m_Volume->GetPlaneCoord( 0, i );
	
	Vector3D v( X, Y, Z );
	v -= centerAbsolute;

	for ( int dim = 0; dim < 3; ++dim )
	  {
	  v[dim] /= radiusAbsolute[dim];
	  }

	if ( v.EuclidNorm() <= 1.0 )
	  this->m_Volume->SetDataAt( value, offset );
	}
      }
    }
}

void
cmtk::UniformVolumePainter::DrawBox
( const IntROI3D& box, const Types::DataItem value )
{
  for ( int k = box.From[2]; k < box.To[2]; ++k )
    {
    for ( int j = box.From[1]; j < box.To[1]; ++j )
      {
      for ( int i = box.From[0]; i < box.To[0]; ++i )
	{
	this->m_Volume->SetDataAt( value, this->m_Volume->GetOffsetFromIndex( i, j, k ) );
	}
      }
    }
}
