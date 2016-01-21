/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 2014 SRI International
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

#include "cmtkPolynomialXform.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

bool
PolynomialXform::ApplyInverse
( const Self::SpaceVectorType& v, Self::SpaceVectorType& u, const Types::Coordinate accuracy ) const
{
  return this->ApplyInverseWithInitial( v, u, v*this->GetGlobalAffineMatrix().GetInverse(), accuracy );
}


const CoordinateMatrix3x3 
PolynomialXform::GetJacobian( const Self::SpaceVectorType& v ) const
{
  const Self::SpaceVectorType& vRel = v - this->m_Center;

  CoordinateMatrix3x3 J = CoordinateMatrix3x3::Identity();

  size_t paramIdx = 0;
  for ( size_t monomialIdx = 0; monomialIdx < this->m_NumberOfMonomials; ++monomialIdx )
    {
    const Types::Coordinate monomialValueDX = Polynomial<4,Types::Coordinate>::EvaluateMonomialDXAt( monomialIdx, vRel[0], vRel[1], vRel[2] );
    const Types::Coordinate monomialValueDY = Polynomial<4,Types::Coordinate>::EvaluateMonomialDYAt( monomialIdx, vRel[0], vRel[1], vRel[2] );
    const Types::Coordinate monomialValueDZ = Polynomial<4,Types::Coordinate>::EvaluateMonomialDZAt( monomialIdx, vRel[0], vRel[1], vRel[2] );

    for ( size_t i = 0; i < 3; ++i, ++paramIdx )
      {
      J[0][i] += this->m_Parameters[paramIdx] * monomialValueDX;
      J[1][i] += this->m_Parameters[paramIdx] * monomialValueDY;
      J[2][i] += this->m_Parameters[paramIdx] * monomialValueDZ;
      }
    }

  return J;
}

const CoordinateMatrix3x3
PolynomialXform::GetLinearMatrix() const
{
  CoordinateMatrix3x3 m3 = CoordinateMatrix3x3::Identity();

  // There is a linear matrix only of degree is non-zero, i.e., the
  // polynomial transform includes at least the linear monomials.
  if (m_Degree > 0) {
    size_t paramIdx = 3;
    for ( size_t j = 0; j < 3; ++j ) {
      for ( size_t i = 0; i < 3; ++i, ++paramIdx ) {
	m3[j][i] += this->m_Parameters[paramIdx];
      }
    }
  }

  return m3;
}

const
AffineXform::MatrixType
PolynomialXform::GetGlobalAffineMatrix() const
{
  // initialize top-left 3x3 as the linear matrix.
  CoordinateMatrix3x3 m3( this->GetLinearMatrix() );

  // transform center using linear matrix
  const Self::SpaceVectorType cM = this->m_Center * m3;
  
  // initialize top-left 3x3 as the linear matrix.
  AffineXform::MatrixType m4x4( this->GetLinearMatrix() );

  // fill in translation parameters, accounting for center.
  for ( size_t i = 0; i < 3; ++i )
    {
    m4x4[3][i] = this->m_Parameters[i] - cM[i] + this->m_Center[i];
    }

  return m4x4;
}

//@}
} // namespace cmtk
