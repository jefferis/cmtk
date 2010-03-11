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

  const int* dims = iMap->GetDims();
#pragma omp parallel for
  for ( int k = 2; k < dims[2]-2; ++k )
    {
    for ( int j = 2; j < dims[1]-2; ++j )
      {
      for ( int i = 2; i < dims[0]-2; ++i )
	{
	Matrix3x3<Types::DataItem> hessian;
	iMap->GetHessianAt( hessian, i, j, k );
	
	EigenSystemSymmetricMatrix3x3<Types::DataItem> eigenSystem( hessian, false /*sortAbsolute*/ );
	
	Types::DataItem result = iMap->GetDataAt( i, j, k );
	if ( eigenSystem.GetNthEigenvalue( 2-this->m_Dimensionality ) < 0 )
	  {
	  Vector3D gradient;
	  iMap->GetGradientAt( gradient, i, j, k );

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
