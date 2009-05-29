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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
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
  v.XYZ[0] -= Origin[0];
  v.XYZ[1] -= Origin[1];
  v.XYZ[2] -= Origin[2];

  p.XYZ[0] = 
    ( v.XYZ[0] * DirectionX[0] + v.XYZ[1] * DirectionX[1] + v.XYZ[2] * DirectionX[2] ) /
    ( MathUtil::Square( DirectionX[0] ) + MathUtil::Square( DirectionX[1] ) + MathUtil::Square( DirectionX[2] ) );
  p.XYZ[1] = 
    ( v.XYZ[0] * DirectionY[0] + v.XYZ[1] * DirectionY[1] + v.XYZ[2] * DirectionY[2] ) /
    ( MathUtil::Square( DirectionY[0] ) + MathUtil::Square( DirectionY[1] ) + MathUtil::Square( DirectionY[2] ) );
  p.XYZ[2] = 0;
}

void
Plane::ProjectPixel( const Vector3D& v, unsigned int& i, unsigned int& j ) const
{
  Vector3D q(v), p;
  this->Project( p, q );

  i = MathUtil::Round( p.XYZ[0] / Spacing[0] );
  j = MathUtil::Round( p.XYZ[1] / Spacing[1] );
}

} // namespace cmtk
