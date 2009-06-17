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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkInfinitePlane.h>

#include <cmtkMathUtil.h>

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

InfinitePlane::InfinitePlane()
{
  Origin.Set( 0, 0, 0 );
  Rho = Theta = Phi = 0;
  this->Update();
}

void
InfinitePlane::Update()
{
  Types::Coordinate radTheta = static_cast<Types::Coordinate>( MathUtil::DegToRad( Theta ) );
  Types::Coordinate radPhi = static_cast<Types::Coordinate>( MathUtil::DegToRad( Phi ) );

  Normal.XYZ[0] = cos( radTheta ) * sin( radPhi );
  Normal.XYZ[1] = sin( radTheta ) * sin( radPhi );
  Normal.XYZ[2] = cos( radPhi );

  SquareNormal = Normal * Normal;
}

void
InfinitePlane::SetNormal( const Vector3D& normal )
{
  this->Normal = (1.0 / normal.EuclidNorm()) * normal;
  
  const Types::Coordinate radPhi = acos( this->Normal.XYZ[2] );
  this->Phi = static_cast<Types::Coordinate>( MathUtil::RadToDeg( radPhi ) );

  const Types::Coordinate sinPhi = sin( radPhi );
  if ( sinPhi != 0 )
    this->Theta = static_cast<Types::Coordinate>( MathUtil::RadToDeg( asin( this->Normal.XYZ[1] / sinPhi ) ) );
  else
    this->Theta = 0;

  this->SquareNormal = this->Normal * this->Normal;
}

AffineXform* 
InfinitePlane::GetAlignmentXform( const byte normalAxis ) const
{
  Types::Coordinate angles[3] = { 0, 0, 0 };
  Types::Coordinate xlate[3] = { 0, 0, 0 };
  
  AffineXform *alignment = new AffineXform;

  switch ( normalAxis ) 
    {
    // YZ plane, i.e., normal to X
    case 0:
    {
    // first make y component zero.
    angles[2] = static_cast<Types::Coordinate>( -MathUtil::RadToDeg( atan2( Normal[1], Normal[0]  ) ) );
    
    // compute x component of normal vector after first rotation; remember that y component will be zero after this rotation.
    const Types::Coordinate newNormal0 = MathUtil::Sign( Normal[0] ) * sqrt( 1 - Normal[2]*Normal[2] );
    angles[1] = static_cast<Types::Coordinate>( -MathUtil::RadToDeg( atan2( Normal[2], newNormal0 ) ) );
    break;
    }

    // XZ plane, normal to Y
    case 1:
    {
    // first make x component zero.
    angles[2] = static_cast<Types::Coordinate>( -MathUtil::RadToDeg( atan2( Normal[0], Normal[1]  ) ) );
    
    // compute y component of normal vector after first rotation; remember that x component will be zero after this rotation.
    const Types::Coordinate newNormal1 = MathUtil::Sign( Normal[1] ) * sqrt( 1 - Normal[2]*Normal[2] );
    angles[0] = static_cast<Types::Coordinate>( -MathUtil::RadToDeg( atan2( Normal[2], newNormal1 ) ) );
    break;
    }

    // XY plane, normal to Z
    case 2:
    {
    // first make x component zero.
    angles[1] = static_cast<Types::Coordinate>( -MathUtil::RadToDeg( atan2( Normal[0], Normal[2]  ) ) );
    
    // compute z component of normal vector after first rotation; remember that x component will be zero after this rotation.
    const Types::Coordinate newNormal2 = MathUtil::Sign( Normal[2] ) * sqrt( 1 - Normal[1]*Normal[1] );
    angles[0] = static_cast<Types::Coordinate>( -MathUtil::RadToDeg( atan2( Normal[1], newNormal2 ) ) );
    break;
    }
    }
  
  alignment->ChangeCenter( this->GetOrigin().XYZ );
  alignment->SetAngles( angles );
  
  xlate[normalAxis] = Rho;
  alignment->SetXlate( xlate );

  return alignment;
}

} // namespace cmtk
