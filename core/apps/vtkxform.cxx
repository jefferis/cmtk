/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
//  $Revision: 1652 $
//
//  $LastChangedDate: 2010-05-14 14:45:52 -0700 (Fri, 14 May 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>

#include <Base/cmtkXform.h>
#include <Base/cmtkXformList.h>

#include <IO/cmtkXformIO.h>
#include <IO/cmtkXformListIO.h>
#include <IO/cmtkVolumeIO.h>

#include <iostream>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <iterator>

int
doMain( const int argc, const char* argv[] )
{
  cmtk::Types::Coordinate inversionTolerance = 0.001;

  std::vector<std::string> inputXformPaths;

  const char* sourceImagePath = NULL;
  const char* targetImagePath = NULL;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Apply coordinate transformation to point coordinates in VTK file." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "An ASCII-format VTK file is read from standard input and a user-provided coordinate transformation (optionally inverted) is applied to the vertex coordinates.  A VTK file with transformed points is then written to standard output." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "vtkxform [options] transformation" );      

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( "inversion-tolerance" ), &inversionTolerance, "Numerical tolerance of B-spline inversion in mm. Smaller values will lead to more accurate inversion, but may increase failure rate." );

    cl.AddOption( Key( "source-image" ), &sourceImagePath, "Set source image of the transformation (i.e., the image that the transformation maps points FROM) to correct for differences in orientation and coordinate space." );
    cl.AddOption( Key( "target-image" ), &targetImagePath, "Set target image of the transformation (i.e., the image that the transformation maps points TO) to correct for differences in orientation and coordinate space." );

    cl.AddParameterVector( &inputXformPaths, "XformList", "List of concatenated transformations. Insert '--inverse' to use the inverse of the transformation listed next." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );  

    cl.Parse( argc, argv );
    }
  catch ( cmtk::CommandLine::Exception ex )
    {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::XformList xformList = cmtk::XformListIO::MakeFromStringList( inputXformPaths );
  xformList.SetEpsilon( inversionTolerance );

  if ( sourceImagePath )
    {
    cmtk::UniformVolume::SmartConstPtr sourceImage( cmtk::VolumeIO::ReadOriented( sourceImagePath ) );
    if ( ! sourceImage )
      {
      cmtk::StdErr << "ERROR: could not read source image '" << sourceImagePath << "'\n";
      throw cmtk::ExitException( 1 );
      }
    xformList.AddToFront( cmtk::AffineXform::SmartPtr( new cmtk::AffineXform( sourceImage->GetImageToPhysicalMatrix() ) )->GetInverse() );
    }
  
  if ( targetImagePath )
    {
    cmtk::UniformVolume::SmartConstPtr targetImage( cmtk::VolumeIO::ReadOriented( targetImagePath ) );
    if ( ! targetImage )
      {
      cmtk::StdErr << "ERROR: could not read target image '" << targetImagePath << "'\n";
      throw cmtk::ExitException( 1 );
      }
    xformList.Add( cmtk::AffineXform::SmartPtr( new cmtk::AffineXform( targetImage->GetImageToPhysicalMatrix() ) ) );
    }
  
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
      
      if ( ! valid )
	{
	// well, not sure what to do now... we should delete the current point from the
	// mesh, but updating the connectivity isn't a local operation. We could also
	// keep track of the previous and the next point and put the failed one in the
	// middle, but that would require a memory and all kinds of special case treatment
	// (multiple consecutive failues, failures at either end, ...) Also assumes
	// a 1D mesh (polyline). So maybe not.
	}

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

#include "cmtkSafeMain"
