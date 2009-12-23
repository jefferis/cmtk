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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkconfig.h>

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>

#include <cmtkXform.h>
#include <cmtkVector3D.h>

#include <cmtkXformIO.h>

#include <cmtkXformList.h>

#include <stdio.h>
#include <list>

bool Help = false;
bool Verbose = false;

std::list<cmtk::Vector3D> LocationList;

void
AddProbeLocation( const char* argv )
{
  double x, y, z;
  if ( 3 == sscanf( argv, "%lf,%lf,%lf", &x, &y, &z ) ) 
    {
    LocationList.push_back( cmtk::Vector3D( x, y, z ) );
    } 
  else
    {
    cmtk::StdErr << "WARNING: '" << argv << "' is not a valid location (must be x,y,z)\n";
    }
}

cmtk::XformList XformList;

int 
main( const int argc, const char* argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Probe transformation at user-defined coordinates" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] [--inverse] xform0 [[--inverse] xform1 ...]" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );
    cl.AddCallback( Key( 'p', "probe" ), AddProbeLocation, "Probe location" );

    cl.Parse();
    
    const char* next = cl.GetNext();
    while ( next ) 
      {
      const bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
      if ( inverse )
	next = cl.GetNext();
      
      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
      if ( ! xform )
	{
	cmtk::StdErr << "ERROR: could not read transformation from " << next << "\n";
	exit( 1 );
	}
      
      XformList.Add( xform, inverse );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception e ) 
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }
  
  for ( std::list<cmtk::Vector3D>::iterator probeIt = LocationList.begin(); probeIt != LocationList.end(); ++probeIt ) 
    {
    cmtk::Vector3D v( *probeIt );
    cmtk::StdOut.printf( "%lf %lf %lf -> ", (double)probeIt->XYZ[0], (double)probeIt->XYZ[1], (double)probeIt->XYZ[2] );
    if ( XformList.ApplyInPlace( v ) )
      cmtk::StdOut.printf( "%lf %lf %lf\n", (double)v.XYZ[0], (double)v.XYZ[1], (double)v.XYZ[2] );
    else
      cmtk::StdOut.printf( "FAILED\n" );
    }
  
  return 0;
}

