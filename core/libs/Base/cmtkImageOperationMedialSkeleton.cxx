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

cmtk::UniformVolume::SmartPtr  
cmtk::ImageOperationMedialSkeleton
::Apply( cmtk::UniformVolume::SmartPtr& volume )
{
  UniformVolume::SmartPtr iMap( new DistanceMapType( volume, DistanceMapType::INSIDE ) );

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
	Types::Coordinate result = iMap->GetDataAt( i, j, k );

	const Types::Coordinate distance = result;
	for ( int kk = -1; (result>0) && (kk < 2); ++kk )
	  for ( int jj = -1; (result>0) && (jj < 2); ++jj )
	    for ( int ii = -1; (result>0) && (ii < 2); ++ii )
	      {
	      if ( kk || jj || ii )
		{
		const Types::Coordinate delta = sqrt( MathUtil::Square( ii * deltas[0] ) + MathUtil::Square( jj * deltas[1] ) +MathUtil::Square( kk * deltas[2] ) );
		if ( (iMap->GetDataAt( i+ii, j+jj, k+kk ) - distance - delta) >= -1e-5 )
		  result = 0;
		}
	      }
	skeleton->SetDataAt( result, i, j, k );
	}
      }
    }
  
  return skeleton;
}
