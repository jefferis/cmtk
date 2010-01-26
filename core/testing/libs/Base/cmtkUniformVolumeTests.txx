/*
//
//  Copyright 2004-2010 SRI International
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkUniformVolume.h>

// test "GridMatches" function
int
testUniformVolumeMatches()
{
  const int dims1[3] = { 10, 11, 12 };
  const float size1a[3] = { 9, 10, 11 };
  const float size1b[3] = { 10, 11, 12 };
  
  const int dims2[3] = { 11, 12, 10 };
  const float size2[3] = { 11, 12, 10 };
  
  cmtk::UniformVolume volume1a( dims1, size1a );
  cmtk::UniformVolume volume1b( dims1, size1b );
  cmtk::UniformVolume volume1c( dims1, size1b );

  cmtk::UniformVolume volume2( dims2, size2 );

  if ( volume1a.GridMatches( volume1b ) )
    return 1;

  if ( !volume1b.GridMatches( volume1c ) )
    return 1;

  if ( volume1a.GridMatches( volume2 ) )
    return 1;
  
  return 0;
}
