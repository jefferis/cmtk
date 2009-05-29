/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkconfig.h>

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>

#include <cmtkXform.h>
#include <cmtkVector3D.h>

#include <cmtkXformIO.h>

#include <stdio.h>
#include <list>

bool Help = false;
bool Verbose = false;

std::list<cmtk::Vector3D> LocationList;

const char* 
AddProbeLocation( const char* argv )
{
  double x, y, z;
  if ( 3 == sscanf( argv, "%lf,%lf,%lf", &x, &y, &z ) ) {
    LocationList.push_back( cmtk::Vector3D( x, y, z ) );
  } else {
    cmtk::StdErr << "WARNING: '" << argv
	      << "' is not a valid location (must be x,y,z)\n";
  }
  return NULL;
}

std::list<cmtk::Xform::SmartPtr> XformList;

int 
main( const int argc, const char* argv[] )
{
  cmtk::CommandLine::PrintProgramID( argv[0], __DATE__, __TIME__ );

  try {
    cmtk::CommandLine cl( argc, argv );

    cl.AddSwitch( Key( 'h', "help" ), &Help, true, "Print command line help" );
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddCallback( Key( 'p', "probe", AddProbeLocation, "Probe location" ) );

    cl.Parse();

    if ( Help ) {
      cl.PrintHelp( "Extended volume reformatter tool.",
		    "[options] xform [xform ...]" );
      exit( 2 );
    }

    const char* next = cl.GetNext();
    while ( next ) {
      bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
      if ( inverse ) next = cl.GetNext();

      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
      if ( ! xform ) {
	cmtk::StdErr << "ERROR: could not read transformation from "
		  << next << "\n";
	exit( 1 );
      }

      XformList.push_back( xform );

      next = cl.GetNextOptional();
    }
  }
  catch ( cmtk::CommandLine::Exception e ) {
    cmtk::StdErr << e << "\n";
    exit( 1 );
  }

  unsigned int xformIdx = 0;

  std::list<cmtk::Xform::SmartPtr>::iterator xformIt = XformList.begin();
  for ( ; xformIt != XformList.end(); ++xformIt, ++xformIdx ) {

    printf( "Transformation #%d:\n", xformIdx );

    std::list<cmtk::Vector3D>::iterator probeIt = LocationList.begin();
    for ( ; probeIt != LocationList.end(); ++probeIt ) {
      cmtk::Vector3D u( *probeIt ), v;
      v = (*xformIt)->Apply( u );
      printf( "\t%f %f %f -> %f %f %f\n", 
	      u.XYZ[0], u.XYZ[1], u.XYZ[2],
	      v.XYZ[0], v.XYZ[1], v.XYZ[2] );
    }

  }
  
  return 0;
}

