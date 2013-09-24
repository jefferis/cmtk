/*
//
//  Copyright 1997-2012 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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
#include <System/cmtkConsole.h>

#include <Base/cmtkXform.h>
#include <Base/cmtkXformList.h>

#include <IO/cmtkXformIO.h>
#include <IO/cmtkXformListIO.h>
#include <IO/cmtkVolumeIO.h>

#include <iostream>
#include <sstream>

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
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Apply coordinate transformations to point coordinates in .fib file." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "A file with fiber tracking results from the UNC Fiber Tracking tool is read from Standard Input and one or more (concatenated) coordinate transformations are applied to all fiber point coordinates. "
		       "The result is written to Standard Output, again in UNC fiber file format." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( "inversion-tolerance" ), &inversionTolerance, "Numerical tolerance of B-spline inversion in mm. Smaller values will lead to more accurate inversion, but may increase failure rate." );

    cl.AddOption( Key( "source-image" ), &sourceImagePath, "Set source image of the transformation (i.e., an image defining the space in which fiber tracking was performed) "
		  "to correct for differences in orientation and coordinate space." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddOption( Key( "target-image" ), &targetImagePath, "Set target image of the transformation (i.e., the image that the fiber track points are mapped into) "
		  "to correct for differences in orientation and coordinate space." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

    cl.AddParameterVector( &inputXformPaths, "XformList", "List of concatenated transformations. Insert '--inverse' to use the inverse of the transformation listed next."
			   "(If the first transformation in the sequence is inverted, then '--inverse' must be preceded by '--', i.e., use '-- --inverse xform.path')." )->SetProperties( cmtk::CommandLine::PROPS_XFORM )
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OPTIONAL );  
    
    cl.Parse( argc, argv );
    }
  catch ( cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::XformList xformList = cmtk::XformListIO::MakeFromStringList( inputXformPaths );
  xformList.SetEpsilon( inversionTolerance );

  if ( sourceImagePath )
    {
    cmtk::UniformVolume::SmartPtr sourceImage( cmtk::VolumeIO::ReadOriented( sourceImagePath ) );
    if ( ! sourceImage )
      {
      cmtk::StdErr << "ERROR: could not read source image '" << sourceImagePath << "'\n";
      throw cmtk::ExitException( 1 );
      }
    
    sourceImage->ChangeCoordinateSpace( sourceImage->GetMetaInfo( cmtk::META_SPACE_ORIGINAL ) );
    xformList.AddToFront( cmtk::AffineXform::SmartPtr( new cmtk::AffineXform(sourceImage->GetImageToPhysicalMatrix() ) ) );
    }

  if ( targetImagePath )
    {
    cmtk::UniformVolume::SmartPtr targetImage( cmtk::VolumeIO::ReadOriented( targetImagePath ) );
    if ( ! targetImage )
      {
      cmtk::StdErr << "ERROR: could not read target image '" << targetImagePath << "'\n";
      throw cmtk::ExitException( 1 );
      }

    targetImage->ChangeCoordinateSpace( targetImage->GetMetaInfo( cmtk::META_SPACE_ORIGINAL ) );
    xformList.Add( cmtk::AffineXform::SmartPtr( new cmtk::AffineXform( targetImage->GetImageToPhysicalMatrix().GetInverse() ) ) );
    }
  
  cmtk::Xform::SpaceVectorType xyz;    
  std::vector<std::string> outputPointLines;

  std::string line;
  while ( !std::cin.eof() )
    {
    std::getline( std::cin, line );

    if ( line.compare( 0, 7, "NPoints" ) != 0 )
      {
      std::cout << line << std::endl;
      }
    else
      {
      const size_t npoints = atoi( line.substr( line.find( '=' )+1, std::string::npos ).c_str() );
      outputPointLines.resize( npoints );

      size_t npointsOut = 0;

      // skip "Points = "
      std::getline( std::cin, line );
      std::cout << line << "\n";
      
      for ( size_t n = 0; n<npoints; ++n )
	{
	// read x,y,z from beginning of line
	std::cin >> xyz[0] >> xyz[1] >> xyz[2];
	
	// read everything else on the same line
	std::string restOfLine;
	std::getline( std::cin, restOfLine );
	
	// Apply transformation sequence
	const bool valid = xformList.ApplyInPlace( xyz );
	
	if ( valid )
	  {
	  std::stringstream sout;
	  sout << xyz[0] << " " << xyz[1] << " " << xyz[2];;
	  outputPointLines[npointsOut++] = std::string( sout.str() + " " + restOfLine );
	  }
	}
      
      std::cout << "NPoints = " << outputPointLines.size() << "\n";
      std::cout << "Points =\n";
      for ( size_t n = 0; n < npointsOut; ++n )
	{
	std::cout << outputPointLines[n] << "\n";
	}
      }
    }

  // if we got here, the program probably ran
  return 0;
}

#include "cmtkSafeMain"
