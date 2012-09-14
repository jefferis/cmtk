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
//  $Revision: 2398 $
//
//  $LastChangedDate: 2010-10-05 14:54:37 -0700 (Tue, 05 Oct 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/
#include <iostream>

#include <Base/cmtkMathUtil.h>

// test whether uniform random numbers are reasonably uniform; does NOT test whether they are random!
int
testMathUtilUniformRandom()
{
  const size_t nSamples = 100000;
  const size_t nBins = 50;

  unsigned int histogram[nBins];
  memset( histogram, 0, sizeof( histogram ) );

  for ( size_t n = 0; n < nSamples; ++n )
    {
    ++histogram[static_cast<size_t>( nBins * cmtk::MathUtil::UniformRandom() )];
    }

  // allow +/- 5% deviation in each bin. This seems to be what Numerical Recipes can do.
  const unsigned int lower = static_cast<unsigned int>( 0.95 * nSamples / nBins );
  const unsigned int upper = static_cast<unsigned int>( 1.05 * nSamples / nBins );

  unsigned int countOutside = 0;
  for ( size_t n = 0; n < nBins; ++n )
    {    
    if ( (histogram[n] < lower ) || (histogram[n] > upper ) )
      ++countOutside;
    }

  const unsigned int threshold = static_cast<unsigned int>( 0.05 * nBins );
  if ( countOutside > threshold )
    {
    std::cerr << "Too many bins outside +/- 5% range. Actual: " << countOutside << ", threshold: " << threshold << std::endl;
    return 1;
    }
  
  return 0;
}
