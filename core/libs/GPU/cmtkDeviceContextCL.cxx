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

#include "DeviceContextCL.h"

#include <System/Exception.h>

cmtk::DeviceContextCL::DeviceContextCL()
{
  cl_int error = CL_SUCCESS;
  this->m_Context = clCreateContextFromType( 0, CL_DEVICE_TYPE_ALL, NULL, NULL, &error );

  if ( error != CL_SUCCESS )
    {
    throw Exception( "clCreateContextFromType() failed" );
    }

  size_t nDevices = 0;
  error = clGetContextInfo( this->m_Context, CL_CONTEXT_DEVICES, 0, NULL, &nDevices );
  if ( error != CL_SUCCESS )
    {
    throw Exception( "clGetContextInfo() failed" );
    }

  this->m_DeviceIDs.resize( nDevices );

  error = clGetContextInfo( this->m_Context, CL_CONTEXT_DEVICES, nDevices, &this->m_DeviceIDs[0], NULL );  
  if ( error != CL_SUCCESS )
    {
    throw Exception( "clGetContextInfo() failed" );
    }
}

cmtk::DeviceContextCL::~DeviceContextCL()
{
  clReleaseContext( this->m_Context );
}
