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

//  const int crop[] = {127,146,17,130,149,20};
//  iMap->SetCropRegion( crop, crop+3 );
//  iMap = UniformVolume::SmartPtr( iMap->GetCroppedVolume() );

  UniformVolume::SmartPtr skeleton( iMap->CloneGrid() );
  skeleton->CreateDataArray( TYPE_COORDINATE );
  skeleton->GetData()->ClearArray();

  const Types::Coordinate* deltas = iMap->GetDelta();

  const int* dims = iMap->GetDims();
#pragma omp parallel for
  for ( int k = 1; k < dims[2]-2; ++k )
    {
    for ( int j = 1; j < dims[1]-2; ++j )
      {
      for ( int i = 1; i < dims[0]-2; ++i )
	{
	Types::DataItem n[3][3][3];
	
	for ( int kk = 0; kk < 3; ++kk )
	  for ( int jj = 0; jj < 3; ++jj )
	    for ( int ii = 0; ii < 3; ++ii )
	      {
	      n[kk][jj][ii] = iMap->GetDataAt( i+ii, j+jj, k+kk );
	      }
	
	Matrix3x3<Types::DataItem> hessian;
	hessian[0][0] = (n[0][0][2] - 2*n[0][0][1] + n[0][0][0] ) / (deltas[0]*deltas[0]);
	hessian[1][2] = (n[0][2][0] - 2*n[0][1][0] + n[0][0][0] ) / (deltas[1]*deltas[1]);
	hessian[1][2] = (n[2][0][0] - 2*n[1][0][0] + n[0][0][0] ) / (deltas[2]*deltas[2]);
	
	hessian[0][1] = hessian[1][0] =  (n[0][1][1] - n[0][0][1] - n[0][1][0] + n[0][0][0] ) / (deltas[0]*deltas[1]);
	hessian[0][2] = hessian[2][0] =  (n[1][0][1] - n[0][0][1] - n[1][0][0] + n[0][0][0] ) / (deltas[0]*deltas[2]);
	hessian[1][2] = hessian[2][1] =  (n[1][1][0] - n[0][1][0] - n[1][0][0] + n[0][0][0] ) / (deltas[1]*deltas[2]);
	
	EigenSystemSymmetricMatrix3x3<Types::DataItem> eigenSystem( hessian, false /*sortAbsolute*/ );
	
	Types::DataItem result = n[1][1][1];
	if ( eigenSystem.GetNthEigenvalue( 2-this->m_Dimensionality ) < 0 )
	  {
	  Vector3D gradient( (n[0][0][1] - n[0][0][0]) / deltas[0], (n[0][1][0] - n[0][0][0]) / deltas[1], (n[1][0][0] - n[0][0][0]) / deltas[2] );
	  for ( int n = 0; n < 3-this->m_Dimensionality; ++n )
	    {
	    Vector3D ev;
	    eigenSystem.GetNthEigenvector( n, ev.XYZ );

	    if ( fabs( ev * gradient ) > 1e-2 )
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
