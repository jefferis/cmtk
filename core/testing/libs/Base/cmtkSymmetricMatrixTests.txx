/*
//
//  Copyright 2010-2011 SRI International
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

#include <Base/cmtkSymmetricMatrix.h>

#include <System/cmtkConsole.h>

// test symmetric behaviour
int
testSymmetricMatrix()
{
  cmtk::SymmetricMatrix<float> m( 2 );

  m(0,0) = -1;
  m(1,1) = -2;

  if ( m(0,0) != -1 )
    {
    cmtk::StdErr << "testSymmetricMatrix test failed at label A\n";
    return 1;
    }

  if ( m(1,1) != -2 )
    {
    cmtk::StdErr << "testSymmetricMatrix test failed at label B\n";
    return 1;
    }

  m(0,1) = 0;
  m(1,0) = 1;
  
  if ( m(0,1) != 1 )
    {
    cmtk::StdErr << "testSymmetricMatrix test failed at label C\n";
    return 1;
    }

  m(0,1) = 2;
  if ( m(1,0) != 2 )
    {
    cmtk::StdErr << "testSymmetricMatrix test failed at label D\n";
    return 1;
    }

  if ( m(0,0) != -1 )
    {
    cmtk::StdErr << "testSymmetricMatrix test failed at label E\n";
    return 1;
    }

  if ( m(1,1) != -2 )
    {
    cmtk::StdErr << "testSymmetricMatrix test failed at label F\n";
    return 1;
    }
  
  return 0;
}

// test resize functionality
int
testSymmetricMatrixResize()
{
  cmtk::SymmetricMatrix<float> m( 1 );

  m(0,0) = 1;
  
  m.Resize( 2 );

  m(0,1) = 2;
  m(1,1) = 3;

  if ( m(0,0) != 1 )
    {
    cmtk::StdErr << "testSymmetricMatrixResize failed at label A\n";
    return 1;
    }

  if ( m(1,0) != 2 )
    {
    cmtk::StdErr << "testSymmetricMatrixResize failed at label B\n";
    return 1;
    }

  if ( m(1,1) != 3 )
    {
    cmtk::StdErr << "testSymmetricMatrixResize failed at label C\n";
    return 1;
    }

  m.Resize( 1 );

  if ( m(0,0) != 1 )
    {
    cmtk::StdErr << "testSymmetricMatrixResize failed at label D\n";
    return 1;
    }

  return 0;
}

// test equality/inequality operators
int
testSymmetricMatrixEqual()
{
  cmtk::SymmetricMatrix<float> m1( 1 );
  cmtk::SymmetricMatrix<float> m2( 2 );

  if ( m1 == m2 )
    {
    cmtk::StdErr << "testSymmetricMatrixEqual failed at label A\n";
    return 1;
    }

  m1.Resize( 2 );
  m1(0,0) = m2(0,0) = 1;
  m1(1,0) = m2(1,0) = 2;
  m1(1,1) = m2(1,1) = 3;

  if ( m1 != m2 )
    {
    cmtk::StdErr << "testSymmetricMatrixEqual failed at label B\n";
    return 1;
    }

  return 0;
}
