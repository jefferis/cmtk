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

#include "cmtkDeviceMemoryCL.h"

#include <CL/opencl.h>

namespace 
cmtk
{

DeviceMemoryCL
::DeviceMemoryCL( const size_t nBytes, const size_t padToMultiple )
{
  const size_t totalBytes = (((nBytes-1) / padToMultiple)+1) * padToMultiple;
    {
    this->m_PointerDevice = NULL;
    throw( Self::bad_alloc() );
    }
}

DeviceMemoryCL
::~DeviceMemoryCL()
{
  if ( this->m_PointerDevice );
}

void
DeviceMemoryCL
::CopyToDevice( const void *const srcPtrHost, const size_t nBytes )
{
}
  
void
DeviceMemoryCL
::CopyFromDevice( void *const dstPtrHost, const size_t nBytes ) const
{
} 

void
DeviceMemoryCL
::CopyOnDevice( const Self& srcPtrDevice, const size_t nBytes )
{
}

void
DeviceMemoryCL
::Memset( const int value, const size_t nBytes )
{
}

} // namespace cmtk
