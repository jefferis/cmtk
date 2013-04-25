/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
  cmtk::Types::Coordinate inversionTolerance = 1e-8;

  std::vector<std::string> inputXformPaths;

  const char* sourceImagePath = NULL;
  const char* targetImagePath = NULL;

  const char* separator = " ";
  int precision = 9;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Apply coordinate transformation to point coordinates from text stream." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "An ASCII-format list of point coordinates is read from standard input and a user-provided sequence of coordinate transformations (each optionally inverted) is applied to them. "
		       "The transformed points are then written to standard output." );

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Xform", "Transformation Options" );
    cl.AddOption( Key( "inversion-tolerance" ), &inversionTolerance, "Numerical tolerance of B-spline inversion in mm. Smaller values will lead to more accurate inversion, but may increase failure rate." );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( "separator" ), &separator, "Field separator - this string or character is used to separate x, y, and z coordinates and z from the remainder of every line." );
    cl.AddOption( Key( "precision" ), &precision, "Floating point precision for output file." );
    cl.EndGroup();

    cl.BeginGroup( "Spaces", "Image Space Options" );    
    cl.AddOption( Key( "source-image" ), &sourceImagePath, "Set source image of the transformation (i.e., the image that the transformation maps points FROM) to correct for differences in orientation and coordinate space." );
    cl.AddOption( Key( "target-image" ), &targetImagePath, "Set target image of the transformation (i.e., the image that the transformation maps points TO) to correct for differences in orientation and coordinate space." );
    cl.EndGroup();

    cl.AddParameterVector( &inputXformPaths, "XformList", "List of concatenated transformations. Insert '--inverse' to use the inverse of the transformation listed next. "
			   "(If the first transformation in the sequence is inverted, then '--inverse' must be preceded by '--', i.e., use '-- --inverse xform.path')." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );  

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

  std::cout << std::setprecision( precision );
  
  cmtk::Xform::SpaceVectorType xyz;
  std::string restOfLine;
  while ( std::cin )
    {
    std::cin >> xyz[0] >> xyz[1] >> xyz[2];
    if ( std::cin )
      {
      std::getline( std::cin, restOfLine );
      
      // Apply transformation sequence
      const bool valid = xformList.ApplyInPlace( xyz );
      
      std::cout << xyz[0] << separator << xyz[1] << separator << xyz[2];
      
      if ( ! valid )
	{
	std::cout << " FAILED";
	}
      
      std::cout << separator << restOfLine << std::endl;
      }
    }
  
  // if we got here, the program probably ran
  return 0;
}

#include "cmtkSafeMain"
