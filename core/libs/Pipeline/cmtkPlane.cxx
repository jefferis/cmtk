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

#include <cmtkconfig.h>

#include <cmtkPlane.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

Plane::Plane()
{
  Dims[0] = Dims[1] = 0;
  Spacing[0] = Spacing[1] = 1;
  Origin[0] = Origin[1] = Origin[2] = 0;
  DirectionX[0] = 1;
  DirectionY[1] = 1;
  DirectionX[1] = DirectionX[2] = DirectionY[0] = DirectionY[2] = 0;
}

void Plane::CopyStructure( const Plane *plane )
{
  this->SetDims( plane->GetDims() );
  this->SetSpacing( plane->GetSpacing() );
  this->SetOrigin( plane->GetOrigin() );
  this->SetDirectionX( plane->GetDirectionX() );
  this->SetDirectionY( plane->GetDirectionY() );
}

void Plane::Project( Vector3D& p, const Vector3D& q ) const
{
  Vector3D v( q );
  v[0] -= Origin[0];
  v[1] -= Origin[1];
  v[2] -= Origin[2];

  p[0] = 
    ( v[0] * DirectionX[0] + v[1] * DirectionX[1] + v[2] * DirectionX[2] ) /
    ( MathUtil::Square( DirectionX[0] ) + MathUtil::Square( DirectionX[1] ) + MathUtil::Square( DirectionX[2] ) );
  p[1] = 
    ( v[0] * DirectionY[0] + v[1] * DirectionY[1] + v[2] * DirectionY[2] ) /
    ( MathUtil::Square( DirectionY[0] ) + MathUtil::Square( DirectionY[1] ) + MathUtil::Square( DirectionY[2] ) );
  p[2] = 0;
}

void
Plane::ProjectPixel( const Vector3D& v, unsigned int& i, unsigned int& j ) const
{
  Vector3D q(v), p;
  this->Project( p, q );

  i = MathUtil::Round( p[0] / Spacing[0] );
  j = MathUtil::Round( p[1] / Spacing[1] );
}

} // namespace cmtk
