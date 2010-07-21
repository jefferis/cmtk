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

#include "cmtkAnatomicalOrientation.h"

#include <string>
#include <iostream>

namespace
cmtk
{

/** \addtogroup Base */
//@{

AnatomicalOrientation::PermutationMatrix::PermutationMatrix
( const FixedVector<3,int>& sourceDims, const std::string& curOrientation, const char newOrientation[3] ) 
{
  // Build a permutation matrix and store it in compressed form 
  for ( int i = 0; i < 3; i++ )
    {
    for ( int j = 0; j < 3; j++ )
      {
      if ( newOrientation[i] == curOrientation[j] )
        {
        this->m_Axes[i] = j; 
        this->m_Multipliers[i] = 1;
        this->m_Offsets[i] = 0;
        break;
        }
      else if ( AnatomicalOrientation::OnSameAxis( newOrientation[i], curOrientation[j] ) )
        {
        this->m_Axes[i] = j; 
        this->m_Multipliers[i] = -1;
        this->m_Offsets[i] = sourceDims[j] - 1;
        break;
        }
      }
    }

  this->m_NewDims = this->GetPermutedArray( sourceDims );
}

AffineXform::MatrixType
AnatomicalOrientation::PermutationMatrix::GetPermutedMatrix( const AffineXform::MatrixType& inMatrix ) const
{
  AffineXform::MatrixType outMatrix, permutation;

  for ( int j = 0; j < 3; ++j )
    {
    for ( int i = 0; i < 3; ++i )
      {
      if ( i == this->m_Axes[j] )
	permutation[i][j] = this->m_Multipliers[j];
      else
	permutation[i][j] = 0;
      }
    permutation[3][j] = this->m_Offsets[j];
    }
  outMatrix = permutation.Invert() * inMatrix;
  
  return outMatrix;
}

} // namespace cmtk
