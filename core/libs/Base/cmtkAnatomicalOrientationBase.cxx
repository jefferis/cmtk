/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkAnatomicalOrientationBase.h"

#include <cassert>
#include <stdlib.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

const char *const AnatomicalOrientationBase::ORIENTATION_STANDARD = "RAS";

const char *const AnatomicalOrientationBase::SPACE_CMTK = "RAS";
const char *const AnatomicalOrientationBase::SPACE_ITK = "LPS";

const char* 
AnatomicalOrientationBase
::GetClosestOrientation( const char* desiredOrientation, const char *const availableOrientations[] )
{
  const char* result = NULL;
  int minPenalty = 100;

  const char *const *next = availableOrientations;
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
AnatomicalOrientationBase::OnSameAxis( const char from, const char to )
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
