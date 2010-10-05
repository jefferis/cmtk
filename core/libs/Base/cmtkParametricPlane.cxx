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

#include "cmtkParametricPlane.h"

#include <Base/cmtkMathUtil.h>

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

ParametricPlane::ParametricPlane()
  : Rho( 0 ),
    Theta( 0 ),
    Phi( 0 )
{
  this->m_Origin = Self::CoordinateVectorType( Self::CoordinateVectorType::Init( 0 ) );
  this->Update();
}

void
ParametricPlane::Update()
{
  Normal[0] = MathUtil::Cos( Theta ) * MathUtil::Sin( Phi );
  Normal[1] = MathUtil::Sin( Theta ) * MathUtil::Sin( Phi );
  Normal[2] = MathUtil::Cos( Phi );
  
  this->SquareNormal = Normal * Normal;
}

void
ParametricPlane::SetNormal( const Self::CoordinateVectorType& normal )
{
  this->Normal = (1.0 / normal.RootSumOfSquares()) * normal;
  
  this->Phi = MathUtil::ArcCos( this->Normal[2] );

  const Types::Coordinate sinPhi = MathUtil::Sin( this->Phi );
  if ( sinPhi != 0 )
    this->Theta = MathUtil::ArcSin( this->Normal[1] / sinPhi );
  else
    this->Theta = Units::Degrees( 0 );

  this->SquareNormal = this->Normal * this->Normal;
}

AffineXform* 
ParametricPlane::GetAlignmentXform( const byte normalAxis ) const
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
    angles[2] = -1.0 * Units::Degrees( MathUtil::ArcTan2( Normal[1], Normal[0] ) ).Value();
    
    // compute x component of normal vector after first rotation; remember that y component will be zero after this rotation.
    const Types::Coordinate newNormal0 = MathUtil::Sign( Normal[0] ) * sqrt( 1 - Normal[2]*Normal[2] );
    angles[1] = -1.0 * Units::Degrees( MathUtil::ArcTan2( Normal[2], newNormal0 ) ).Value();
    break;
    }

    // XZ plane, normal to Y
    case 1:
    {
    // first make x component zero.
    angles[2] = -1.0 * Units::Degrees( MathUtil::ArcTan2( Normal[0], Normal[1] ) ).Value();
    
    // compute y component of normal vector after first rotation; remember that x component will be zero after this rotation.
    const Types::Coordinate newNormal1 = MathUtil::Sign( Normal[1] ) * sqrt( 1 - Normal[2]*Normal[2] );
    angles[0] = -1.0 * Units::Degrees( MathUtil::ArcTan2( Normal[2], newNormal1 ) ).Value();
    break;
    }

    // XY plane, normal to Z
    case 2:
    {
    // first make x component zero.
    angles[1] = -1.0 * Units::Degrees( MathUtil::ArcTan2( Normal[0], Normal[2] ) ).Value();
    
    // compute z component of normal vector after first rotation; remember that x component will be zero after this rotation.
    const Types::Coordinate newNormal2 = MathUtil::Sign( Normal[2] ) * sqrt( 1 - Normal[1]*Normal[1] );
    angles[0] = -1.0 * Units::Degrees( MathUtil::ArcTan2( Normal[1], newNormal2 ) ).Value();
    break;
    }
    }
  
  alignment->ChangeCenter( this->GetOrigin() );
  alignment->SetAngles( angles );
  
  xlate[normalAxis] = Rho;
  alignment->SetXlate( xlate );

  return alignment;
}

AffineXform::MatrixType
ParametricPlane::GetMirrorXformMatrix() const
{
  // put together zero-offset mirror matrix
  AffineXform::MatrixType M = AffineXform::MatrixType::IdentityMatrix;

  for ( int i = 0; i < 3; ++i ) 
    {
    for ( int j = 0; j < 3; ++j ) 
      {
      M[i][j] -= 2.0 * this->Normal[i]*this->Normal[j] / this->SquareNormal;
      }
    }

  FixedVector<3,Types::Coordinate> mo = this->m_Origin;
  mo *= M;

  for ( int j = 0; j < 3; ++j ) 
    {
    M[3][j] = this->m_Origin[j] - mo[j] + 2 * this->Rho * this->Normal[j] / this->SquareNormal;
    }

  return M;
}

} // namespace cmtk
