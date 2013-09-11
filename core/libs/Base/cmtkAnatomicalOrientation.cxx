/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010, 2013 SRI International
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

#include <System/cmtkConsole.h>

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
AnatomicalOrientation
::GetOrientationFromDirections( char* orientation,  const AffineXform::MatrixType& directions, const char* spaceAxes )
{
  const Types::Coordinate spacing[3] =
    {
      sqrt( directions[0][0]*directions[0][0] + directions[0][1]*directions[0][1] + directions[0][2]*directions[0][2] ),
      sqrt( directions[1][0]*directions[1][0] + directions[1][1]*directions[1][1] + directions[1][2]*directions[1][2] ),
      sqrt( directions[2][0]*directions[2][0] + directions[2][1]*directions[2][1] + directions[2][2]*directions[2][2] )
    };

  // keep track of which axes are already used in the direction code - need 4 entries to allow for access when axes are non-ortogonal
  bool axisUsed[4] = { false, false, false, true };

  for ( int axis = 0; axis < 3; ++axis )
    {
    // skip axes already used
    int maxDim = 0;
    while ( axisUsed[maxDim] ) ++maxDim;
    
    // get closest aligned of remaining axes
    Types::Coordinate max = fabs( directions[axis][0] / spacing[axis] );
    for ( int dim = 1; dim < 3; ++dim )
      {
      const Types::Coordinate positive = fabs( directions[axis][dim] / spacing[axis] );
      if ( (positive > max) && !axisUsed[dim] )
	{
	max = positive;
	maxDim = dim;
	}
      else
	{
	if ( positive == max )
	  {
	  maxDim = 3;
	  }
	}
      }
    
    if ( maxDim == 3 )
      {
      StdErr << "WARNING: image seems to have an oblique orientation. This is not going to end well...\n";
      }
    
    orientation[axis] = (directions[axis][maxDim] > 0) ? spaceAxes[maxDim] : Self::OppositeDirection( spaceAxes[maxDim] );
    axisUsed[maxDim] = true;
    }
  orientation[3] = 0;
}

void
AnatomicalOrientation
::GetImageToSpaceAxesPermutation( int (&imageToSpaceAxesPermutation)[3][3], const char* orientation, const char* spaceAxes )
{
  for ( int j = 0; j < 3; ++j )
    {
    for ( int i = 0; i < 3; ++i )
      {
      if ( orientation[j] == spaceAxes[i] )
	 imageToSpaceAxesPermutation[j][i] = 1;
      else 
	if ( Self::OnSameAxis( orientation[j], spaceAxes[i] ) )
	  imageToSpaceAxesPermutation[j][i] = -1;
	else
	  imageToSpaceAxesPermutation[j][i] = 0;
      }
    }
}

} // namespace cmtk
