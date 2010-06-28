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

// test "DeviceHistogram" class
int
testDeviceHistogramEntropy()
{
  try
    {
    cmtk::DeviceHistogram::SmartPtr histogram = cmtk::DeviceHistogram::Create( 100 );

    float floatHost[100];
    // compute entropy for all-zeros
    for ( size_t i = 0; i < 100; ++i )
      floatHost[i] = 1;
    histogram->GetDataOnDevice().CopyToDevice( floatHost, 100 );
    std::cerr << "Entropy: " << histogram->GetEntropy() << std::endl;

    // compute entropy for 50x 0, 50x 1
    for ( size_t i = 0; i < 50; ++i )
      floatHost[i] = 0;
    histogram->GetDataOnDevice().CopyToDevice( floatHost, 100 );
    std::cerr << "Entropy: " << histogram->GetEntropy() << std::endl;

    // compute entropy for 50x "0 1" alternating
    for ( size_t i = 0; i < 50; ++i )
      {
      floatHost[i<<1] = 1;
      floatHost[1+(i<<1)] = 0;
      }
    histogram->GetDataOnDevice().CopyToDevice( floatHost, 100 );
    std::cerr << "Entropy: " << histogram->GetEntropy() << std::endl;
    }
  catch ( std::bad_alloc )
    {
    std::cerr << "Caught bad_alloc()" << std::endl;
    return 1;
    }

  return 0;
}

