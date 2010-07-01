/*
//
//  Copyright 2004-2010 SRI International
//
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

#include <iostream>
#include <cstring>

#include <cmtkTestFunctionMap.h>

#include "cmtkDeviceHistogramTests.txx"
#include "cmtkDeviceMemoryTests.txx"
#include "cmtkDeviceUniformVolumeTests.txx"

int
main( const int argc, const char* argv[] )
{
  cmtk::TestFunctionMap map;
  map.AddTest( "DeviceHistogramEntropy", &testDeviceHistogramEntropy );
  map.AddTest( "DeviceHistogramPopulate", &testDeviceHistogramPopulate );
  map.AddTest( "DeviceMemory", &testDeviceMemory );
  map.AddTest( "DeviceUniformVolume", &testDeviceUniformVolume );

  // is test name given on command line?
  if ( argc < 2 )
    {
    }
  else
    {
    return map.RunTestByName( argv[1] );
    }
}
