/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2009 SRI International
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

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>

#include <IO/cmtkVolumeIO.h>
#include <Base/cmtkMathUtil.h>

unsigned int NumberOfPoints = 0;
bool UseImageAsMask = false;
bool WriteIndex = false;

const char* ImageFileName = NULL;

int
doMain( const int argc, const char *argv[] )
{
  try
    {
    cmtk::CommandLine cl;
    typedef cmtk::CommandLine::Key Key;
    
    cl.AddSwitch( Key( 'i', "write-index" ), &WriteIndex, true, "For each rancom pixel, write its index, not its location in world coordinates." );
    cl.AddOption( Key( 'n', "num-points" ), &NumberOfPoints, "Number of points to generate." );
    cl.AddSwitch( Key( 'm', "mask" ), &UseImageAsMask, true, "Use image data as a mask: only generates points with non-zero pixels." );

    cl.Parse( argc, argv );

    ImageFileName = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::Read( ImageFileName ) );
  for ( unsigned int i  = 0; i < NumberOfPoints; ++i )
    {
    int x, y, z;
    for ( bool okay = false; !okay; )
      {
      x = static_cast<int>( cmtk::MathUtil::UniformRandom() * volume->GetDims()[0] );
      y = static_cast<int>( cmtk::MathUtil::UniformRandom() * volume->GetDims()[1] );
      z = static_cast<int>( cmtk::MathUtil::UniformRandom() * volume->GetDims()[2] );

      okay = !UseImageAsMask || volume->GetDataAt( x, y, z );
      }
    if ( WriteIndex )
      {
      std::cout << x << "\t" << y << "\t" << z << std::endl;
      }
    else
      {
      std::cout << volume->GetPlaneCoord(cmtk::AXIS_X,x) << "\t" << volume->GetPlaneCoord(cmtk::AXIS_Y,y) << "\t" << volume->GetPlaneCoord(cmtk::AXIS_Z,z) << std::endl;
      }
    }
}

#include "cmtkSafeMain"
