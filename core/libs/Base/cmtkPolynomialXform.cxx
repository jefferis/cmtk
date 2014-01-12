/*
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
      J[i][0] += this->m_Parameters[paramIdx] * monomialValueDX;
      J[i][1] += this->m_Parameters[paramIdx] * monomialValueDY;
      J[i][2] += this->m_Parameters[paramIdx] * monomialValueDZ;
      }
    }

  return J;
}

//@}
} // namespace cmtk
