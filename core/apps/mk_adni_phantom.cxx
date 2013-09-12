/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#include <fstream>

int
doMain( const int argc, const char* argv[] )
{
  cmtk::Types::Coordinate resolution = 1.0;
  bool labels = false;

  std::string outputFileName = "phantom.nii";
  std::string outputLabelsName;
  std::string outputLandmarksName;

  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Generate ADNI phantom image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Generate image of the ADNI structural imaging calibration phantom (a.k.a. Magphan EMR051)." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( "resolution" ), &resolution, "Set output image resolution in [mm]" );
    cl.AddSwitch( Key( "labels" ), &labels, true, "Draw each marker sphere with a label value defined by its index in the marker table. Otherwise, estimated T1 is used." );
    cl.AddOption( Key( "write-labels" ), &outputLabelsName, "Optional path to write text file with label names." )->SetProperties( cmtk::CommandLine::PROPS_OUTPUT );
    cl.AddOption( Key( "write-landmarks" ), &outputLandmarksName, "Optional path to write text file with landmark locations." )->SetProperties( cmtk::CommandLine::PROPS_OUTPUT );

    cl.AddParameter( &outputFileName, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e;
    return 1;
    }

  cmtk::VolumeIO::Write( *(cmtk::MagphanEMR051::GetPhantomImage( resolution, labels )), outputFileName );

  // write optional labels file
  if ( !outputLabelsName.empty() )
    {
    std::ofstream stream( outputLabelsName.c_str() );
    if ( stream.good() )
      {
      for ( size_t i = 0; i < cmtk::MagphanEMR051::NumberOfSpheres; ++i )
	{
	stream << i+1 << "\t" << cmtk::MagphanEMR051::SphereName( i ) << std::endl;
	}
      }
    }
  
  // write optional landmarks file
  if ( !outputLandmarksName.empty() )
    {
    std::ofstream stream( outputLandmarksName.c_str() );
    if ( stream.good() )
      {
      for ( size_t i = 0; i < cmtk::MagphanEMR051::NumberOfSpheres; ++i )
	{
	stream << cmtk::MagphanEMR051::SphereCenter( i )[0] << "\t" << cmtk::MagphanEMR051::SphereCenter( i )[1] << "\t" << cmtk::MagphanEMR051::SphereCenter( i )[2] << "\t" 
	       << cmtk::MagphanEMR051::SphereName( i ) << "\t" << std::endl;
	}
      }
    }
  
  return 0;
}

#include "cmtkSafeMain"
