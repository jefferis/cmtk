/*
//
//  Copyright 2004-2010 SRI International
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

#include <cmtkPGM.h>
#include <cmtkScalarImage.h>

// test whether  we can read 8bit and 16bit PGM files.
int
testReadPGM()
{
  cmtk::ScalarImage::SmartPtr image8( cmtk::PGM::Read( CMTK_DATADIR "/axial.pgm" ) );
  if ( ! image8 )
    {
    std::cerr << "ERROR: could not read 8bit PGM test image 'axial.pgm'" << std::endl;
    return 1;
    }
  
  cmtk::ScalarImage::SmartPtr image16( cmtk::PGM::Read( CMTK_DATADIR "/axial16.pgm" ) );
  if ( ! image16 )
    {
    std::cerr << "ERROR: could not read 16bit PGM test image 'axial16.pgm'" << std::endl;
    return 1;
    }
  
  return 0;
}
