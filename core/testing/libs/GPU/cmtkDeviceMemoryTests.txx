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

#include <cmtkDeviceMemoryCUDA.h>

#include <cuda_runtime_api.h>

// test "DeviceMemoryCUDA" class
int
testDeviceMemoryCUDA()
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
    cmtk::DeviceMemoryCUDA<float>::SmartPtr floatCUDA = cmtk::DeviceMemoryCUDA<float>::Create( 100 );
    floatCUDA->CopyToDevice( floatHost, 100 );
    floatCUDA->CopyFromDevice( floatHost, 100 );

    int intHost[100];
    cmtk::DeviceMemoryCUDA<int>::SmartPtr intCUDA = cmtk::DeviceMemoryCUDA<int>::Create( 100 );
    intCUDA->CopyToDevice( intHost, 100 );
    intCUDA->CopyFromDevice( intHost, 100 );

    char charHost[100];
    cmtk::DeviceMemoryCUDA<char>::SmartPtr charCUDA = cmtk::DeviceMemoryCUDA<char>::Create( 100 );
    charCUDA->CopyToDevice( charHost, 100 );
    charCUDA->CopyFromDevice( charHost, 100 );

    cmtk::DeviceMemoryCUDA<float>::SmartPtr float2CUDA = cmtk::DeviceMemoryCUDA<float>::Create( 100 );
    float2CUDA->CopyOnDevice( *floatCUDA, 100 );

    size_t memFreeAfter, memTotalAfter;
    if ( cudaMemGetInfo( &memFreeAfter, &memTotalAfter ) != cudaSuccess )
      {
      std::cerr << "Second call to cudaMemGetInfo() failed." << std::endl;
      return 1;
      }

    if ( memFreeBefore == memFreeAfter )
      {
      std::cerr << "Free CUDA device memory constant despite allocation at" << memFreeBefore << std::endl;
      return 1;
      }  
    }
  catch ( cmtk::DeviceMemoryBaseCUDA::bad_alloc )
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
    std::cerr << "Total CUDA device memory changed by " << (memTotalBefore - memTotalAfter) << std::endl;
    return 1;
    }
  
  if ( memFreeBefore != memFreeAfter )
    {
    std::cerr << "Free CUDA device memory changed by" << (memFreeBefore - memFreeAfter) << std::endl;
    return 1;
    }
  
  return 0;
}

