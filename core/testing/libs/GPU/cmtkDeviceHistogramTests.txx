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

#include <cmtkDeviceHistogram.h>

#include <cuda_runtime_api.h>
#include <cmath>

int
checkEntropy( const std::string& testName, const float hData[100], cmtk::DeviceHistogram& dHist, const float baseline )
{
  dHist.GetDataOnDevice().CopyToDevice( hData, 100 );
  const float entropy = dHist.GetEntropy();

  if ( fabs( entropy - baseline ) > 1e-5 )
    {
    std::cerr << "Test " << testName << " entropy " << entropy << " deviates from baseline " << baseline << std::endl;
    return 1;
    }
  return 0;
}

// test "DeviceHistogram" class
int
testDeviceHistogramEntropy()
{
  try
    {
    cmtk::DeviceHistogram::SmartPtr histogram100 = cmtk::DeviceHistogram::Create( 100 );
    cmtk::DeviceHistogram::SmartPtr histogram200 = cmtk::DeviceHistogram::Create( 200 );
    
    float floatHost[100];

    // compute entropy for all-zeros
    for ( size_t i = 0; i < 100; ++i )
      floatHost[i] = 0;
    if ( checkEntropy( "AllZeros100", floatHost, *histogram100, 0 ) || checkEntropy( "AllZeros200", floatHost, *histogram200, 0 ) )
      return 1;
    
    // compute entropy for single non-zero bin
    floatHost[0] = 1;
    if ( checkEntropy( "SingleBin100", floatHost, *histogram100, 0 ) || checkEntropy( "SingleBin200", floatHost, *histogram200, 0 ) )
      return 1;

    // compute entropy for all-ones
    for ( size_t i = 0; i < 100; ++i )
      floatHost[i] = 1;
    if ( checkEntropy( "AllOnes100", floatHost, *histogram100, 4.60517 ) || checkEntropy( "AllOnes200", floatHost, *histogram200, 4.60517 ) )
      return 1;

    // compute entropy for 50x 0, 50x 1
    for ( size_t i = 0; i < 50; ++i )
      floatHost[i] = 0;
    if ( checkEntropy( "50One50Zero100", floatHost, *histogram100, 3.91202 ) || checkEntropy( "50One50Zero200", floatHost, *histogram200, 3.91202 ) )
      return 1;

    // compute entropy for 50x "0 1" alternating
    for ( size_t i = 0; i < 50; ++i )
      {
      floatHost[i<<1] = 1;
      floatHost[1+(i<<1)] = 0;
      }
    if ( checkEntropy( "50OneZeroPairs100", floatHost, *histogram100, 3.91202 ) || checkEntropy( "50OneZeroPairs200", floatHost, *histogram200, 3.91202 ) )
      return 1;
    }
  catch ( std::bad_alloc )
    {
    std::cerr << "Caught bad_alloc()" << std::endl;
    return 1;
    }

  return 0;
}

