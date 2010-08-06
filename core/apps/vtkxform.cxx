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
//  $Revision: 1652 $
//
//  $LastChangedDate: 2010-05-14 14:45:52 -0700 (Fri, 14 May 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include "System/cmtkCommandLine.h"
#include "System/cmtkConsole.h"

#include "IO/cmtkXformIO.h"
#include "Registration/cmtkXformList.h"

#include "Base/cmtkXform.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iterator>

int
main( const int argc, const char* argv[] )
{
  bool verbose = false;

  cmtk::Types::Coordinate inversionTolerance = 0.001;
  cmtk::XformList xformList;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Apply coordinate transformation to point coordinates in VTK file (standard input) and write equivalent file with transformed points to standard output." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] transformation" );      

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Print each point to STDERR (as well as stdout)" );
    cl.AddOption( Key( "inversion-tolerance" ), &inversionTolerance, "Numerical tolerance of B-spline inversion in mm. Smaller values will lead to more accurate inversion, but may increase failure rate." );

    cl.Parse( argc, argv );

    const char* next = cl.GetNextOptional();
    while (next)
      {
      bool inverse = false;
      if ( !strcmp( next, "-i" ) || !strcmp( next, "--inverse" ) )
	{
	inverse = true;
	next = cl.GetNext();
	}

      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, verbose ) );
      xformList.Add( xform, inverse );
      next = cl.GetNextOptional();
      }
    }
  catch ( cmtk::CommandLine::Exception ex )
    {
    cmtk::StdErr << ex << "\n";
    exit( 1 );
    }

  xformList.SetEpsilon( inversionTolerance );
  
  // Is VTK file stored in binary format?
  bool binaryMode = false;

  // First, read everything up to and including the "POINTS" line and write everything to output unchanged
  std::string line;
  while ( !std::cin.eof() )
    {
    std::getline( std::cin, line );
    std::cout << line << std::endl;
    
    if ( ! line.compare( 0, 6, "BINARY" ) )
      binaryMode = true;

    if ( ! line.compare( 0, 6, "POINTS" ) )
      break;
    }
  
  // If we're not at EOF, then "line" must be "POINTS <npoints> ..."
  if ( ! std::cin.eof() )
    {
    // Parse number of points out of line
    std::stringstream sstream( line.substr( 7 ) );
    size_t npoints;
    sstream >> npoints;
    
    // Repeat npoints times
    cmtk::Xform::SpaceVectorType xyz;
    for ( size_t n = 0; n<npoints; ++n )
      {
      // Read original point coordinates from file
      if ( binaryMode )
	{
	float xyzFloat[3];
	std::cin.read( reinterpret_cast<char*>( &xyzFloat[0] ), sizeof( xyzFloat ) );

#ifndef WORDS_BIGENDIAN
	for ( size_t i = 0; i<3; ++i )
	  cmtk::Memory::ByteSwapInPlace( xyzFloat[i] );
#endif // #ifndef WORDS_BIGENDIAN

	xyz = cmtk::FixedVector<3,float>( xyzFloat );
	}
      else
	{
	std::cin >> xyz[0] >> xyz[1] >> xyz[2];
	}
      
      // Apply transformation sequence
      const bool valid = xformList.ApplyInPlace( xyz );
      
      // Write transformed point to output
      // Read original point coordinates from file
      if ( binaryMode )
	{
	float xyzFloat[3] = { xyz[0], xyz[1], xyz[2] };

#ifndef WORDS_BIGENDIAN
	for ( size_t i = 0; i<3; ++i )
	  cmtk::Memory::ByteSwapInPlace( xyzFloat[i] );
#endif // #ifndef WORDS_BIGENDIAN

	std::cout.write( reinterpret_cast<const char*>( &xyzFloat[0] ), sizeof( xyzFloat ) );
	}
      else
	{
	std::cout << xyz[0] << " " << xyz[1] << " " << xyz[2] << std::endl;
	}
      }
    }

  // Everything else remains unchanged, so copy from input to output.
  char c = std::cin.get();
  while ( std::cin.good() )
    {
    std::cout << c;
    c = std::cin.get();
    }

  // if we got here, the program probably ran
  return 0;
}

