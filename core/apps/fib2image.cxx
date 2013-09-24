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
  std::string trackingImagePath;
  std::string outputImagePath = "fib2image.nii";

  cmtk::Types::DataItem value = 1;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Draw image from point coordinates of fiber tracks from .fib file." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Fiber tracking results from the UNC Fiber Tracking tool are read from Standard Input and all fiber points are drawn into a 3D image. The result is written in one of the supported image file formats." );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( "value" ), &value, "Value used for drawing fiber points." );
    cl.AddOption( Key( 'o', "output" ), &outputImagePath, "Output image file name." );

    cl.AddParameter( &trackingImagePath, "TrackingImage", "Image defining the grid and space in which fiber tracking was performed to correct for differences in orientation and coordinate space." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    
    cl.Parse( argc, argv );
    }
  catch ( cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }

  // read target image in NATIVE orientation because that's what UNC Fiber Tracking uses
  cmtk::UniformVolume::SmartPtr trackingImage( cmtk::VolumeIO::ReadOriented( trackingImagePath ) );
  if ( ! trackingImage )
    {
    cmtk::StdErr << "ERROR: could not read tracking image '" << trackingImagePath << "'\n";
    throw cmtk::ExitException( 1 );
    }
  
  // don't need data
  trackingImage->SetData( cmtk::TypedArray::SmartPtr( NULL ) );
  cmtk::UniformVolume::SmartPtr outputImage( trackingImage->CloneGrid() );

  trackingImage->ChangeCoordinateSpace( trackingImage->GetMetaInfo( cmtk::META_SPACE_ORIGINAL ) );
  cmtk::AffineXform fromTrackingSpace( trackingImage->GetImageToPhysicalMatrix() );

  // don't need data - make new array
  outputImage->CreateDataArray( cmtk::TYPE_SHORT, true /*setToZero*/ );
  
  cmtk::Xform::SpaceVectorType xyz;
  cmtk::DataGrid::IndexType ijk;

  std::string line;
  while ( !std::cin.eof() )
    {
    std::getline( std::cin, line );

    if ( line.compare( 0, 7, "NPoints" ) == 0 )
      {
      const size_t npoints = atoi( line.substr( line.find( '=' )+1, std::string::npos ).c_str() );

      // skip "Points = "
      std::getline( std::cin, line );
      
      for ( size_t n = 0; n<npoints; ++n )
	{
	// read x,y,z from beginning of line
	std::cin >> xyz[0] >> xyz[1] >> xyz[2];
	
	// read everything else on the same line to skip.
	std::string restOfLine;
	std::getline( std::cin, restOfLine );

	// transform from fib space into image space
	if ( outputImage->GetClosestGridPointIndex( fromTrackingSpace.Apply( xyz ), ijk ) )
	  {
	  outputImage->SetDataAt( value, trackingImage->GetOffsetFromIndex( ijk ) );
	  }
	}
      }
    }
  
  // write output image
  cmtk::VolumeIO::Write( *outputImage, outputImagePath );

  // if we got here, the program probably ran
  return 0;
}

#include "cmtkSafeMain"
