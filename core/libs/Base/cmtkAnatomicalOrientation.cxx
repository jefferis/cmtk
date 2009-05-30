/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkAnatomicalOrientation.h>

#include <cmtkConsole.h>

#include <math.h>
#include <assert.h>

#include <string>
#include <sstream>

namespace
cmtk
{

/** \addtogroup Base */
//@{

const char* AnatomicalOrientation::ORIENTATION_STANDARD = "RAS";

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

  for ( int axis = 0; axis < 3; ++axis )
    {
    Types::Coordinate max = fabs( directions[axis][0] / spacing[axis] );
    int maxDim = 0;
    for ( int dim = 1; dim < 3; ++dim )
      {
      const Types::Coordinate positive = fabs( directions[axis][dim] / spacing[axis] );
      if ( positive > max )
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
	if ( OnSameAxis( orientation[j], spaceAxes[i] ) )
	  imageToSpaceAxesPermutation[j][i] = -1;
	else
	  imageToSpaceAxesPermutation[j][i] = 0;
      }
    }
}

const char* 
AnatomicalOrientation
::GetClosestOrientation( const char* desiredOrientation, const char* availableOrientations[] )
{
  const char* result = NULL;
  int minPenalty = 100;

  const char** next = availableOrientations;
  while ( *next )
    {
    int penalty = 0;
    for ( int axis = 0; axis < 3; ++axis )
      {
      if ( desiredOrientation[axis] != (*next)[axis] )
	{
	if ( Self::OnSameAxis( desiredOrientation[axis], (*next)[axis] ) )
	  penalty += 1;
	else
	  penalty += 4;
	}
      }

    if ( penalty < minPenalty )
      {
      result = *next;
      minPenalty = penalty;
      }

    ++next;
    }
  return result;
}

bool
AnatomicalOrientation::OnSameAxis( const char from, const char to )
{
  // Set up lists such that the direction corresponding to the 
  // character at index "from-'A'" is the character that represents
  // the same axis in reverse direction. Lowercase characters are
  // for padding and orientation.

  assert( (from=='A') || (from=='P') || (from=='L') || (from=='R') || (from=='I') || (from=='S') );
  assert( (to=='A') || (to=='P') || (to=='L') || (to=='R') || (to=='I') || (to=='S') );

  return (Self::OppositeDirection( from ) == to);
}

} // namespace cmtk
