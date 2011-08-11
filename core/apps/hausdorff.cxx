/*
//
//  Copyright 2011 SRI International
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

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <IO/cmtkVolumeIO.h>

#include <Registration/cmtkHausdorffDistance.h>

int
doMain
( const int argc, const char *argv[] )
{
  const char* imagePath0 = NULL;
  const char* imagePath1 = NULL;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Hausdorff Distance." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool computes the Hausdorff distance between two label images." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddParameter( &imagePath0, "Image0", "First image path." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &imagePath1, "Image1", "Second image path." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartConstPtr image0 = cmtk::VolumeIO::Read( imagePath0 );
  if ( ! image0 || ! image0->GetData() )
    {
    cmtk::StdErr << "ERROR: unable to read image " << imagePath0 << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartConstPtr image1 = cmtk::VolumeIO::Read( imagePath1 );
  if ( ! image1 || ! image1->GetData() )
    {
    cmtk::StdErr << "ERROR: unable to read image " << imagePath1 << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::StdOut << cmtk::HausdorffDistance( image0, image1 ).GetBinary() << "\n";

  return 0;
}

#include "cmtkSafeMain"
