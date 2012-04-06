/*
//
//  Copyright 2012 SRI International
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

#include "cmtkFitRigidToLandmarks.h"

cmtk::FitRigidToLandmarks::FitRigidToLandmarks( const LandmarkPairList& landmarkPairs )
{
  // first, get the centroids in "from" and "to" space
  cmtk::FixedVector<3,cmtk::Types::Coordinate> cFrom( cmtk::FixedVector<3,cmtk::Types::Coordinate>::Init( 0 ) );
  cmtk::FixedVector<3,cmtk::Types::Coordinate> cTo( cmtk::FixedVector<3,cmtk::Types::Coordinate>::Init( 0 ) );
  
  size_t nLandmarks = 0;
  for ( LandmarkPairList::const_iterator it = landmarkPairs.begin(); it != landmarkPairs.end(); ++it )
    {
    cFrom += it->m_Location;
    cTo += it->m_TargetLocation;
    ++nLandmarks;
    }
  
  cFrom /= nLandmarks;
  cTo /= nLandmarks;
  
  // now compute the transformation matrix for rotation, scale, shear, using the previously computed centroids for reference
  Matrix2D<double> U( 3, 3 ); 
  U.SetAllToZero();
  
  // build the two 3x3 matrices of (t*xT)(x*xT) on the fly.
  for ( LandmarkPairList::const_iterator it = landmarkPairs.begin(); it != landmarkPairs.end(); ++it )
    {
    const cmtk::FixedVector<3,cmtk::Types::Coordinate> x = it->m_Location - cFrom;
    const cmtk::FixedVector<3,cmtk::Types::Coordinate> t = it->m_TargetLocation - cTo;
    
    for ( size_t j = 0; j < 3; ++j )
      {
      for ( size_t i = 0; i < 3; ++i )
	{
	U[i][j] += t[j] * x[i];
	}
      }
    }  

  Matrix2D<double> V( 3, 3 );
  std::vector<double> W( 3 );
  MathUtil::SVD( U, 3, 3, W, V );

  Matrix3x3<Types::Coordinate> matrix;
  for ( size_t j = 0; j < 3; ++j )
    {
    for ( size_t i = 0; i < 3; ++i )
      {
      matrix[j][i] = 0;
      for ( size_t k = 0; k < 3; ++k )
	{
	matrix[j][i] += V[i][k] * U[j][k];
	}
      }
    }

  // put everything together
  AffineXform::MatrixType matrix4x4( matrix );
  this->m_RigidXform = AffineXform::SmartPtr( new AffineXform( matrix4x4 ) );
  this->m_RigidXform->SetTranslation( (cTo - cFrom) );
  this->m_RigidXform->SetCenter( cFrom );
}
