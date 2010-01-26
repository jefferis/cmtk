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

#include <cmtkDataGrid.h>

// test "GridMatches" function
int
testDataGridMatches()
{
  const int dims1[3] = { 10, 11, 12 };
  const int dims2[3] = { 11, 12, 10 };
  
  cmtk::DataGrid grid1a( dims1 );
  cmtk::DataGrid grid1b( dims1 );
  cmtk::DataGrid grid2( dims2 );

  if ( !grid1a.GridMatches( grid1b ) )
    return 1;

  if ( !grid1b.GridMatches( grid1a ) )
    return 1;

  if ( grid1a.GridMatches( grid2 ) )
    return 1;
  
  return 0;
}
