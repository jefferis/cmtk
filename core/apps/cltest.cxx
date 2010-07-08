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

#include <cmtkconfig.h>

#include <CL/opencl.h>

#include <cstdlib>
#include <iostream>
#include <vector>

int
main( const int argc, const char*[] )
{
  cl_int error = CL_SUCCESS;

  cl_context_properties props[] = { CL_CONTEXT_PLATFORM, NULL };
  cl_context context = clCreateContextFromType( props, CL_DEVICE_TYPE_ALL, NULL, NULL, &error );

  if ( error != CL_SUCCESS )
    {
    std::cerr << "clCreateContextFromType() failed with error " << error << std::endl;
    exit( 1 );
    }

  size_t nDevices = 0;
  error = clGetContextInfo( context, CL_CONTEXT_DEVICES, 0, NULL, &nDevices );
  if ( error != CL_SUCCESS )
    {
    std::cerr << "clGetContextInfo() failed" << std::endl;
    exit( 1 );
    }
  
  std::vector<cl_device_id> deviceIDs( nDevices );
  
  error = clGetContextInfo( context, CL_CONTEXT_DEVICES, nDevices, &deviceIDs[0], NULL );  
  if ( error != CL_SUCCESS )
    {
    std::cerr << "clGetContextInfo() failed" << std::endl;
    exit( 1 );
    }

  for ( size_t id = 0; id < deviceIDs.size(); ++id )
    {
    std::cerr << std::endl << "Device #" << id << ":" << std::endl;

    size_t result;
    error = clGetDeviceInfo( deviceIDs[id], CL_DEVICE_MAX_COMPUTE_UNITS, NULL, NULL, &result );
    std::cerr << "\tMax compute units count: " << result << std::endl;
    }

#if 0
  for ( int device = 0; device < deviceCount; ++device )
    {
    std::cerr << std::endl << "Device #" << device << ":" << std::endl;
    struct cudaDeviceProp props;
    if ( cudaGetDeviceProperties( &props, device ) != cudaSuccess )
      {
      std::cerr << "\tFailed to get device properties." << std::endl;
      }
    else
      {
      std::cerr << "\tName: " << props.name << std::endl << std::endl;
      std::cerr << "\tMultiprocessor count: " << props.multiProcessorCount << std::endl;
      std::cerr << "\tCompute capability: " << props.major << "." << props.minor << std::endl;
      std::cerr << "\tTotal memory: " << props.totalGlobalMem << std::endl;
      std::cerr << "\tConstant memory: " << props.totalConstMem << std::endl;
      std::cerr << "\tShared memory per block: " << props.sharedMemPerBlock << std::endl << std::endl;
      std::cerr << "\tWarp size: " << props.warpSize << std::endl;
      std::cerr << "\tMax threads per block: " << props.maxThreadsPerBlock << std::endl;
      std::cerr << "\tMaximum thread block size: (" << props.maxThreadsDim[0] << "," << props.maxThreadsDim[1] << "," << props.maxThreadsDim[2] << ")" << std::endl;
      std::cerr << "\tMaximum grid size size: (" << props.maxGridSize[0] << "," << props.maxGridSize[1] << "," << props.maxGridSize[2] << ")" << std::endl;
      }
    }
#endif

  // if we got here, the program probably ran
  return 0;
}

