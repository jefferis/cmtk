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

#include <Base/cmtkUniformVolume.h>
#include <IO/cmtkVolumeIO.h>

int
main( const int argc, const char*[] )
{
  if ( argc > 3 )
    {
    const cmtk::UniformVolume::IndexType::ValueType dims[3] = {2,2,2};
    const cmtk::UniformVolume::CoordinateVectorType::ValueType size[3] = {1,1,1};
    cmtk::UniformVolume::SmartPtr volume( new cmtk::UniformVolume( cmtk::UniformVolume::IndexType( dims ), cmtk::UniformVolume::CoordinateVectorType( size ) ) );
    }

  // if we got here, the program probably ran
  return 0;
}

