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

#include <Numerics/ap.h>
#include <Numerics/sevd.h>

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

  const Types::Coordinate oneOverN = 1.0 / landmarkPairs.size();
  for ( size_t j = 0; j < 3; ++j )
    {
    for ( size_t i = 0; i < 3; ++i )
      {
      U[i][j] *= oneOverN;
      }
    }

  const Types::Coordinate traceU = U[0][0] + U[1][1] + U[2][2];
  const Types::Coordinate delta[3] = { U[1][2] - U[2][1], U[2][0] - U[0][2], U[0][1] - U[1][0] };

  ap::real_2d_array a;
  a.setbounds( 0, 3, 0, 3 );
  a(0,0) = traceU;

  a(0,1) = a(1,0) = delta[0];
  a(0,2) = a(2,0) = delta[1];
  a(0,3) = a(3,0) = delta[2];

  for ( size_t j = 0; j < 3; ++j )
    {
    for ( size_t i = 0; i < 3; ++i )
      {
      a(1+i,1+j)=a(1+j,1+i) = U[i][j]+U[j][i];
      }
    a(1+j,1+j) -= traceU;
    }
  
  ap::real_1d_array eVals;
  eVals.setbounds( 0, 3 );

  ap::real_2d_array eVecs;
  eVecs.setbounds( 0, 3, 0, 3 );

  if ( smatrixevd( a, 4, true /*ZNeeded*/, true /*isUpper*/, eVals, eVecs ) )
    {
    // eVals/eVecs are sorted in ascending order; use last
    const Types::Coordinate q[4] = { eVecs(3,0),  eVecs(3,1),  eVecs(3,2),  eVecs(3,3) };
    
    // put everything together
    AffineXform::MatrixType matrix;
    matrix[0][0] = q[0]*q[0]+q[1]*q[1]-q[2]*q[2]-q[3]*q[3];
    matrix[0][1] = 2*(q[1]*q[2]-q[0]*q[3]);
    matrix[0][2] = 2*(q[1]*q[3]+q[0]*q[2]);
    matrix[1][0] = 2*(q[1]*q[2]+q[0]*q[3]);
    matrix[1][1] = q[0]*q[0]+q[2]*q[2]-q[1]*q[1]-q[3]*q[3];
    matrix[1][2] = 2*(q[2]*q[3]-q[0]*q[1]);
    matrix[2][0] = 2*(q[1]*q[3]-q[0]*q[2]);
    matrix[2][1] = 2*(q[2]*q[3]+q[0]*q[1]);
    matrix[2][2] = q[0]*q[0]+q[3]*q[3]-q[1]*q[1]-q[2]*q[2];
    matrix[3][3] = 1;
    matrix[3][0] = matrix[3][1] = matrix[3][2] = 0;
    matrix[0][3] = matrix[0][3] = matrix[0][3] = 0;

    this->m_RigidXform = AffineXform::SmartPtr( new AffineXform( matrix ) );
    const AffineXform::SpaceVectorType cFromR = this->m_RigidXform->Apply( cFrom );
    this->m_RigidXform->SetTranslation( (cTo - cFromR) );
    this->m_RigidXform->SetCenter( cFrom );
    }
}
