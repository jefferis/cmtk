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

#include <cmtkMatrix3x3.h>
#include <cmtkEigenSystemSymmetricMatrix3x3.h>

#include <iostream>
#include <math.h>

template<class T>
bool
testMatrixEigensystem
( const T m[3][3], const T evals[3], const T evecs[3][3] )
{
  const T tolerance = 1e-6;
  
  cmtk::Matrix3x3<T> matrix( &m[0][0] );
  cmtk::EigenSystemSymmetricMatrix3x3<T> es( matrix );
  
  // compare eigenvalues
  for ( size_t i = 0; i<3; ++i )
    {
    if ( fabs( evals[i] - es.GetNthEigenvalue(i) ) > tolerance )
      {
      std::cerr << "Eigenvalues do not match." << std::endl
		<< "  ACTUAL: " << es.GetNthEigenvalue(0) << " " << es.GetNthEigenvalue(1) << " " << es.GetNthEigenvalue(2) << std::endl
		<< "  BASELN: " << evals[0] << " " << evals[1] << " " << evals[2] << std::endl;
      return false;
      }
    }

  // compare eigenvectors
  for ( size_t i = 0; i<3; ++i )
    {
    const cmtk::FixedVector<3,T> actual = es.GetNthEigenvector( i );
    for ( size_t j = 0; j<3; ++j )
      {
      if ( fabs( evecs[i][j] - actual[j] ) > tolerance )
	{
	std::cerr << "Eigenvectors do not match." << std::endl;
	for ( size_t ii = 0; ii<3; ++ii )
	  {
	  const cmtk::FixedVector<3,T> ev = es.GetNthEigenvector( ii );
	  std::cerr << "  ACTUAL: " << ev[0] << " " << ev[1] << " " << ev[2] << std::endl;
	  std::cerr << "  BASELN: " << evecs[ii][0] << " " << evecs[ii][1] << " " << evecs[ii][2] << " " << std::endl;
	  }
	return false;
	}
      }
    }
  
  return true;
}

// test eigenvalues computation
int
testEigenSystemSymmetricMatrix3x3()
{
  const double dm1[3][3] = { { 1, 0, 0 }, { 0, -2, 0 }, { 0, 0, 3 } };
  const double evals1[3] = { 1, -2, 3 };
  const double evecs1[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
  if ( ! testMatrixEigensystem( dm1, evals1, evecs1 ) )
    {
    return 1;
    }
  
  return 0;
}
