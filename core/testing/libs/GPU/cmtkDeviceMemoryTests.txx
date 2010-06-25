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

#include <cmtkDeviceMemory.h>

#include <cuda_runtime_api.h>

// test "DeviceMemory" class
int
testDeviceMemory()
{
  size_t memFreeBefore, memTotalBefore;
  if ( cudaMemGetInfo( &memFreeBefore, &memTotalBefore ) != cudaSuccess )
    {
    std::cerr << "Call to cudaMemGetInfo() failed." << std::endl;
    return 1;
    }

  try
    {
    float floatHost[100];
    cmtk::DeviceMemory<float>::SmartPtr floatDevice = cmtk::DeviceMemory<float>::Create( 100 );
    floatDevice->CopyToDevice( floatHost, 100 );
    floatDevice->CopyFromDevice( floatHost, 100 );

    int intHost[100];
    cmtk::DeviceMemory<int>::SmartPtr intDevice = cmtk::DeviceMemory<int>::Create( 100 );
    intDevice->CopyToDevice( intHost, 100 );
    intDevice->CopyFromDevice( intHost, 100 );

    char charHost[100];
    cmtk::DeviceMemory<char>::SmartPtr charDevice = cmtk::DeviceMemory<char>::Create( 100 );
    charDevice->CopyToDevice( charHost, 100 );
    charDevice->CopyFromDevice( charHost, 100 );

    cmtk::DeviceMemory<float>::SmartPtr float2Device = cmtk::DeviceMemory<float>::Create( 100 );
    float2Device->CopyOnDevice( *floatDevice, 100 );

    size_t memFreeAfter, memTotalAfter;
    if ( cudaMemGetInfo( &memFreeAfter, &memTotalAfter ) != cudaSuccess )
      {
      std::cerr << "Second call to cudaMemGetInfo() failed." << std::endl;
      return 1;
      }

    if ( memFreeBefore == memFreeAfter )
      {
      std::cerr << "Free device memory constant despite allocation at" << memFreeBefore << std::endl;
      return 1;
      }  
    }
  catch ( std::bad_alloc )
    {
    std::cerr << "Caught bad_alloc()" << std::endl;
    return 1;
    }

  size_t memFreeAfter, memTotalAfter;
  if ( cudaMemGetInfo( &memFreeAfter, &memTotalAfter ) != cudaSuccess )
    {
    std::cerr << "Third call to cudaMemGetInfo() failed." << std::endl;
    return 1;
    }

  if ( memTotalBefore != memTotalAfter )
    {
    std::cerr << "Total device memory changed by " << (memTotalBefore - memTotalAfter) << std::endl;
    return 1;
    }
  
  if ( memFreeBefore != memFreeAfter )
    {
    std::cerr << "Free device memory changed by" << (memFreeBefore - memFreeAfter) << std::endl;
    return 1;
    }
  
  return 0;
}

