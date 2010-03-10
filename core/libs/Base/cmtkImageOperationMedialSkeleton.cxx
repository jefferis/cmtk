/*
//
//  Copyright 2009,2010 SRI International
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

#include <cmtkImageOperationMedialSkeleton.h>

#include <cmtkEigenSystemSymmetricMatrix3x3.h>

cmtk::UniformVolume::SmartPtr  
cmtk::ImageOperationMedialSkeleton
::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  UniformVolume::SmartPtr iMap( new DistanceMapType( *volume, DistanceMapType::INSIDE ) );

  UniformVolume::SmartPtr skeleton( volume->CloneGrid() );
  skeleton->CreateDataArray( TYPE_COORDINATE );
  skeleton->GetData()->ClearArray();

  const Types::Coordinate* deltas = volume->GetDelta();

  const int* dims = volume->GetDims();
#pragma omp parallel for
  for ( int k = 1; k < dims[2]-1; ++k )
    {
    for ( int j = 1; j < dims[1]-1; ++j )
      {
      for ( int i = 1; i < dims[0]-1; ++i )
	{
	Types::DataItem n[3][3][3];
	
	for ( int kk = -1; kk < 2; ++kk )
	  for ( int jj = -1; jj < 2; ++jj )
	    for ( int ii = -1; ii < 2; ++ii )
	      {
	      n[1+kk][1+jj][1+ii] = iMap->GetDataAt( i+ii, j+jj, k+kk );
	      }
	
	Types::DataItem dx[3][3][2];
	Types::DataItem dy[3][3][2];
	Types::DataItem dz[3][3][2];
	for ( int kk = 0; kk < 3; ++kk )
	  for ( int jj = 0; jj < 3; ++jj )
	    for ( int ii = 0; ii < 2; ++ii )
	      {
	      dx[kk][jj][ii] = (n[kk][jj][ii+1] - n[kk][jj][ii] ) / deltas[0];
	      dy[kk][jj][ii] = (n[kk][ii+1][jj] - n[kk][ii][jj] ) / deltas[1];
	      dz[kk][jj][ii] = (n[ii+1][kk][jj] - n[ii][kk][jj] ) / deltas[2];
	      }
	
	Matrix3x3<Types::DataItem> hessian;
	hessian[0][0] = (dx[1][1][1] - dx[1][1][0]) / deltas[0];
	hessian[1][1] = (dy[1][1][1] - dy[1][1][0]) / deltas[1];
	hessian[2][2] = (dz[1][1][1] - dz[1][1][0]) / deltas[2];
	
	hessian[0][1] = hessian[1][0] = (dx[1][0][1] + dx[1][0][0] + dx[1][2][1] + dx[1][2][0] - 2 * ( dx[1][1][1] + dx[1][1][0] ) ) / (2*deltas[1]);
	hessian[0][2] = hessian[2][0] = (dx[0][1][1] + dx[0][1][0] + dx[2][1][1] + dx[2][1][0] - 2 * ( dx[1][1][1] + dx[1][1][0] ) ) / (2*deltas[2]);
	hessian[1][2] = hessian[2][1] = (dy[0][1][1] + dy[0][1][0] + dy[2][1][1] + dy[2][1][0] - 2 * ( dy[1][1][1] + dy[1][1][0] ) ) / (2*deltas[2]);
	
	EigenSystemSymmetricMatrix3x3<Types::DataItem> eigenSystem( hessian, false /*sortAbsolute*/ );
	
	Types::DataItem result = n[1][1][1];
	if ( eigenSystem.GetNthEigenvalue( 2 - this->m_Dimensionality ) < 0 )
	  {
	  Vector3D gradient( (dx[1][1][0] + dx[1][1][1]) / 2, (dy[1][1][0] + dy[1][1][1]) / 2, (dz[1][1][0] + dz[1][1][1]) / 2 );
	  for ( int n = 0; n < 2 - this->m_Dimensionality; ++n )
	    {
	    Vector3D ev;
	    eigenSystem.GetNthEigenvector( n, ev.XYZ );

	    if ( fabs( ev * gradient ) > 1e-5 )
	      result = 0;
	    }
	  }
	else
	  {
	  result = 0;
	  }
	skeleton->SetDataAt( result, i, j, k );
	}
      }
    }
  
  return skeleton;
}
