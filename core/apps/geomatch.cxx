/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
//  $Revision: 4818 $
//
//  $LastChangedDate: 2013-09-10 11:28:54 -0700 (Tue, 10 Sep 2013) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkFixedSquareMatrix.h>

#include <IO/cmtkVolumeIO.h>

#include <stdio.h>
#include <vector>
#include <string>

cmtk::UniformVolume::SmartPtr
readVolume( const std::string& path, const char* readOrientation )
{
  cmtk::UniformVolume::SmartPtr volume;
  try
    {
    if ( readOrientation )
      volume = cmtk::VolumeIO::ReadOriented( path, readOrientation );
    else
      volume = cmtk::VolumeIO::Read( path );
    }
  catch (...) 
    {
    if ( readOrientation )
      volume = cmtk::VolumeIO::ReadGridOriented( path, readOrientation );
    else
      volume = cmtk::VolumeIO::ReadGrid( path );
    }
  
  if ( ! volume )
    {
    cmtk::StdErr << "ERROR: cannot read image from " << path << "\n";
    throw cmtk::ExitException( 1 );
    }

  return volume;
}

int
doMain( int argc, const char *argv[] )
{
  const char* readOrientation = NULL;
  std::vector<std::string> imagePaths;

  double tolerance = 1e-5;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Check whether the geometries (e.g., grid dimensions, pixel sizes, spatial coordinates) or two or more images match." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool reads two or more images and tests whether their grid dimensions, pixel sizes, and image-to-space transformations match. "
		       "Optionally, all images are reoriented into standard orientation before performing the test. If all images match, the tool returns with exit code 0, otherwise it returns with exit code 2. "
		       "In case of an error (e.g., one of the images can not be read), the exit code is 1." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( "read-ras" ), &readOrientation, "RAS", "Read all images in RAS orientation" );
    cl.AddOption( Key( "tolerance" ), &tolerance, "Numerical tolerance for floating point comparisons (e.g., transformation matrices)" );

    cl.AddParameterVector( &imagePaths, "ImagePaths", "List of image files." );
    
    cl.Parse( argc, const_cast<const char**>( argv ) );
    }
  catch ( const cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }

  // Check if we have enough parameters
  if ( imagePaths.size() < 2 )
    {
    cmtk::StdErr << "ERROR: need at least two image paths.\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartConstPtr firstVolume = readVolume( imagePaths[0], readOrientation );
  for ( size_t i = 1; i < imagePaths.size(); ++i )
    {
    cmtk::UniformVolume::SmartConstPtr nextVolume = readVolume( imagePaths[i], readOrientation );
    if ( ! firstVolume->GridMatches( *nextVolume ) )
      {
      return 2;
      }

    if ( (firstVolume->GetImageToPhysicalMatrix() - nextVolume->GetImageToPhysicalMatrix()).FrobeniusNorm() > tolerance )
      {
      return 2;
      }
    }

  return 0;
}

#include "cmtkSafeMain"

