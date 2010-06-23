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

#include <cmtkUniformVolumeCUDA.h>

// test "UniformVolumeCUDA" class
int
testUniformVolumeCUDA()
{
  try
    {
    const int dims[3] = { 10, 10, 10 };
    const float size[3] = { 9, 9, 9 };

    cmtk::UniformVolume volume( (cmtk::FixedVector<3,int>( dims )), cmtk::FixedVector<3,float>( size ) );

    // first, try to create representation of actually empty volume.
    cmtk::UniformVolumeCUDA::SmartPtr volumeCUDA = cmtk::UniformVolumeCUDA::Create( volume );

    // second, allocate pixel data and create another device instance.
    volume.CreateDataArray( cmtk::TYPE_INT );
    volumeCUDA = cmtk::UniformVolumeCUDA::Create( volume );    

    // third, change pixel data to double precision float and create another device instance.
    volume.CreateDataArray( cmtk::TYPE_DOUBLE );
    volumeCUDA = cmtk::UniformVolumeCUDA::Create( volume );    
    }
  catch ( cmtk::DeviceMemoryBaseCUDA::bad_alloc )
    {
    std::cerr << "Caught bad_alloc()" << std::endl;
    return 1;
    }
  
  return 0;
}

