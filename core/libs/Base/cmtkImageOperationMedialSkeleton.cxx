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

  const Types::Coordinate* delta = volume->GetDelta();

  size_t ofs = 0;
  const int* dims = volume->GetDims();
  for ( int k = 0; k < dims[2]; ++k )
    {
    for ( int j = 0; j < dims[1]; ++j )
      {
      for ( int i = 0; i < dims[0]; ++i, ++ofs )
	{
	const Types::Coordinate distance = iMap->GetDataAt( ofs );
	if ( distance <= 0 )
	  {
	  skeleton->SetDataAt( 0, ofs );
	  }
	else
	  {
	  }
	}
      }
    }
  
  return skeleton;
}
