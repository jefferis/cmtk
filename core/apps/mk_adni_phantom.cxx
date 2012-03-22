/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkMagphanEMR051.h>

#include <IO/cmtkVolumeIO.h>

int
doMain( const int argc, const char* argv[] )
{
  cmtk::Types::Coordinate resolution = 1.0;
  const char* outputFileName = "phantom.nii";

  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Generate ADNI phantom image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Generate image of the ADNI structural imaging calibration phantom (a.k.a. Magphan EMR051)." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( "resolution" ), &resolution, "Set output image resolution in [mm]" );

    cl.AddParameter( &outputFileName, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e;
    return 1;
    }

  cmtk::VolumeIO::Write( *(cmtk::Phantoms::GetPhantomImage(resolution)), outputFileName );
  
  return 0;
}

#include "cmtkSafeMain"
