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

#include "cmtkDeviceTransformedVolumeAxes.h"

#include "System/cmtkException.h"

cmtk::DeviceTransformedVolumeAxes
::DeviceTransformedVolumeAxes( const TransformedVolumeAxes& axesOnHost )
{
  for ( int dim = 0; dim < 3; ++dim )
    {
    const size_t dims = axesOnHost.Dims()[dim];

    std::vector<float> points( 3 * dims );
    for ( size_t idx = 0; idx < dims; ++idx )
      {
      for ( size_t comp = 0; comp < 3; ++comp )
	points[3*idx+comp] = static_cast<float>( axesOnHost[dim][idx][comp] );
      }
    
    this->m_OnDeviceData[dim] = DeviceMemory<float>::Create( 3 * dims );
    this->m_OnDeviceData[dim]->CopyToDevice( &points[0], 3*dims );
    }
}
