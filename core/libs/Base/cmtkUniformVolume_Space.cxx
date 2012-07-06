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

#include "cmtkUniformVolume.h"

#include <Base/cmtkAnatomicalOrientation.h>
#include <Base/cmtkAffineXform.h>

namespace
cmtk
{

void
UniformVolume::CreateDefaultIndexToPhysicalMatrix()
{
  this->m_IndexToPhysicalMatrix = AffineXform::MatrixType::Identity();
  for ( int axis = 0; axis < 3; ++axis )
    for ( int i = 0; i < 3; ++i )
      this->m_IndexToPhysicalMatrix[axis][i] *= this->m_Delta[axis];
}

const UniformVolume::SmartPtr
UniformVolume::GetReoriented( const char* newOrientation ) const
{
  const std::string curOrientation = this->GetMetaInfo( META_IMAGE_ORIENTATION );
  DataGrid::SmartPtr temp( DataGrid::GetReoriented( newOrientation ) );

  AnatomicalOrientation::PermutationMatrix pmatrix( this->m_Dims, curOrientation, newOrientation );
  FixedVector<3,Types::Coordinate> newSize = pmatrix.GetPermutedArray( this->m_Size );
  
  UniformVolume::SmartPtr result( new UniformVolume( temp->GetDims(), newSize, temp->GetData() ) );
  result->m_Offset = pmatrix.GetPermutedArray( this->m_Offset );
  result->m_IndexToPhysicalMatrix = pmatrix.GetPermutedMatrix( this->m_IndexToPhysicalMatrix );

  result->CopyMetaInfo( *temp );
  return result;
}

void
UniformVolume
::ChangeCoordinateSpace( const std::string& newSpace )
{
  const std::string currentSpace = this->GetMetaInfo( META_SPACE );
  if ( currentSpace == newSpace )
    return; // nothing to do.

  int axesPermutation[3][3];
  AnatomicalOrientation::GetImageToSpaceAxesPermutation( axesPermutation, newSpace.c_str(), currentSpace.c_str() );

  AffineXform::MatrixType newMatrix = AffineXform::MatrixType::Identity();
  for ( int j = 0; j < 3; ++j )
    {
    for ( int j2 = 0; j2 < 3; ++j2 )
      {
      if ( axesPermutation[j][j2] != 0 )
	{
	for ( int i = 0; i < 4; ++i )
	  {
	  newMatrix[i][j] = axesPermutation[j][j2] * this->m_IndexToPhysicalMatrix[i][j2];
	  }
	}
      }
    }
  
  this->SetMetaInfo( META_SPACE, newSpace );
  this->m_IndexToPhysicalMatrix = newMatrix;
}

std::string
UniformVolume
::GetOrientationFromDirections() const
{
  const AffineXform::MatrixType& matrix = this->m_IndexToPhysicalMatrix;
  char orientationString[4] = { 0,0,0,0 };
  AnatomicalOrientation::GetOrientationFromDirections( orientationString, matrix, this->GetMetaInfo( META_SPACE ).c_str() );
  return std::string( orientationString );
}

AffineXform::MatrixType
UniformVolume::GetImageToPhysicalMatrix() const
{
  AffineXform::MatrixType matrix = this->m_IndexToPhysicalMatrix;
// mDelta[3] is implicitly == 1 (homogeneous coordinates), so 4th matrix row (translation/coordinate origin) stays untouched
  for ( int i = 0; i < 3; ++i )
    {
    if ( this->m_Delta[i] > 0 )
      {
      for ( int j = 0; j < 3; ++j )
	matrix[i][j] /= this->m_Delta[i];
      }
    }

  return matrix;
}

void
UniformVolume::SetImageToPhysicalMatrix( const AffineXform::MatrixType& matrix )
{
  this->m_IndexToPhysicalMatrix = matrix;
// mDelta[3] is implicitly == 1 (homogeneous coordinates), so 4th matrix row (translation/coordinate origin) stays untouched
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      this->m_IndexToPhysicalMatrix[i][j] *= this->m_Delta[i];
}

} // namespace cmtk
