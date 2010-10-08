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

#include "cmtkUniformVolume.h"

#include <Base/cmtkAnatomicalOrientation.h>
#include <Base/cmtkAffineXform.h>

namespace
cmtk
{

void
UniformVolume::CreateDefaultIndexToPhysicalMatrix()
{
  this->m_IndexToPhysicalMatrix = AffineXform::MatrixType::IdentityMatrix;
  for ( int axis = 0; axis < 3; ++axis )
    for ( int i = 0; i < 3; ++i )
      this->m_IndexToPhysicalMatrix[axis][i] *= this->m_Delta[axis];
}

const UniformVolume::SmartPtr
UniformVolume::GetReoriented( const char* newOrientation ) const
{
  const std::string curOrientation = this->m_MetaInformation[META_IMAGE_ORIENTATION];
  DataGrid::SmartPtr temp( DataGrid::GetReoriented( newOrientation ) );

  AnatomicalOrientation::PermutationMatrix pmatrix( this->m_Dims, curOrientation, newOrientation );
  FixedVector<3,Types::Coordinate> newSize = pmatrix.GetPermutedArray( this->Size );
  
  UniformVolume::SmartPtr result( new UniformVolume( temp->GetDims(), newSize, temp->GetData() ) );
  result->m_Offset = pmatrix.GetPermutedArray( this->m_Offset );
  result->m_IndexToPhysicalMatrix = pmatrix.GetPermutedMatrix( this->m_IndexToPhysicalMatrix, this->Size );

  result->m_MetaInformation = temp->m_MetaInformation;
  return result;
}

void
UniformVolume
::ChangeCoordinateSpace( const std::string& newSpace )
{
  const std::string currentSpace = this->m_MetaInformation[META_SPACE];
  if ( currentSpace == newSpace )
    return; // nothing to do.

  int axesPermutation[3][3];
  AnatomicalOrientation::GetImageToSpaceAxesPermutation( axesPermutation, newSpace.c_str(), currentSpace.c_str() );

  AffineXform::MatrixType newMatrix;
  for ( int j = 0; j < 3; ++j )
    {
    for ( int j2 = 0; j2 < 3; ++j2 )
      {
      if ( axesPermutation[j][j2] )
	{
	for ( int i = 0; i < 4; ++i )
	  {
	  newMatrix[i][j] = axesPermutation[j][j2] * this->m_IndexToPhysicalMatrix[i][j2];
	  }
	}
      }
    }
  
  this->m_MetaInformation[META_SPACE] = newSpace;
  this->m_IndexToPhysicalMatrix = newMatrix;
}

std::string
UniformVolume
::GetOrientationFromDirections() const
{
  const AffineXform::MatrixType& matrix = this->m_IndexToPhysicalMatrix;
  char orientationString[4] = { 0,0,0,0 };
  AnatomicalOrientation::GetOrientationFromDirections( orientationString, matrix, this->m_MetaInformation[META_SPACE].c_str() );
  return std::string( orientationString );
}

AffineXform::MatrixType
UniformVolume::GetImageToPhysicalMatrix() const
{
  AffineXform::MatrixType matrix = this->m_IndexToPhysicalMatrix;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      matrix[i][j] /= this->m_Delta[i];

  return matrix;
}

} // namespace cmtk
