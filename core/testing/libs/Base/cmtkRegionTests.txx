/*
//
//  Copyright 2010 SRI International
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

#include <cmtkFixedVector.h>
#include <cmtkRegion.h>

// test "Size" member function for int regions
int
testRegionSizeInt()
{
  const int regionFromArray[] = { 1, 2, 3 };
  const int regionToArray[] = { 2, 4, 7 };

  const cmtk::Region<3,int> r3( (cmtk::FixedVector<3,int>( regionFromArray )), cmtk::FixedVector<3,int>( regionToArray ) );
  if ( r3.Size() != 8 )
    return 1;

  const cmtk::Region<2,int> r2( (cmtk::FixedVector<2,int>( regionFromArray )), cmtk::FixedVector<2,int>( regionToArray ) );
  if ( r2.Size() != 2 )
    return 1;

  const cmtk::Region<1,int> r1( (cmtk::FixedVector<1,int>( regionFromArray )), cmtk::FixedVector<1,int>( regionToArray ) );
  if ( r1.Size() != 1 )
    return 1;
  
  return 0;
}

// test "Size" member function for float regions
int
testRegionSizeFloat()
{
  const float regionFromArray[] = { 1, 2, 3 };
  const float regionToArray[] = { 2, 4, 7 };

  const cmtk::Region<3,float> r3( (cmtk::FixedVector<3,float>( regionFromArray )), cmtk::FixedVector<3,float>( regionToArray ) );
  if ( r3.Size() != 8 )
    return 1;

  const cmtk::Region<2,float> r2( (cmtk::FixedVector<2,float>( regionFromArray )), cmtk::FixedVector<2,float>( regionToArray ) );
  if ( r2.Size() != 2 )
    return 1;

  const cmtk::Region<1,float> r1( (cmtk::FixedVector<1,float>( regionFromArray )), cmtk::FixedVector<1,float>( regionToArray ) );
  if ( r1.Size() != 1 )
    return 1;
  
  return 0;
}
