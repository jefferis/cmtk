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
#include <Base/cmtkVector3D.h>

#include <IO/cmtkVolumeIO.h>

#include <iostream>

typedef enum
{
  COORDINATES_INDEXED,
  COORDINATES_ABSOLUTE,
  COORDINATES_RELATIVE,
  COORDINATES_PHYSICAL
} CoordinateModeEnum;

int
doMain( const int argc, const char *argv[] )
{
  std::string inputImagePath;
  const char* readOrientation = "RAS";

  CoordinateModeEnum inputMode = COORDINATES_ABSOLUTE;
  CoordinateModeEnum outputMode = COORDINATES_ABSOLUTE;

  const char* radiusStr = "1";

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Search image neighborhoods for pixels." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool reads an image file, as well as a list of pixel coordinates from standard input. For each pixel, a local neighbourhood  in the image is searched for the maximum "
		       "value. The location of the maximum is then written to standard output." );

    typedef cmtk::CommandLine::Key Key;
    cmtk::CommandLine::EnumGroup<CoordinateModeEnum>::SmartPtr inputModeGroup = cl.AddEnum( "input-coordinates", &inputMode, "Coordinate specification mode for program input." );
    inputModeGroup->AddSwitch( Key( "absolute" ), COORDINATES_ABSOLUTE, "Use absolute volume coordinates. For each dimension, the valid range is [0,FOV]." );
    inputModeGroup->AddSwitch( Key( "indexed" ), COORDINATES_INDEXED, "Use grid indexes to specify coordinates. For each dimension, the valid value range is [0,Dims-1]." );
    inputModeGroup->AddSwitch( Key( "relative" ), COORDINATES_RELATIVE, "Use relative volume coordinates. For each dimension, the valid range is [0,1]." );
    inputModeGroup->AddSwitch( Key( "physical" ), COORDINATES_PHYSICAL, "Use physical volume coordinates. "
			       "Each given location is transformed into image coordinates via the inverse of the images's index-to-physical space matrix." );

    cmtk::CommandLine::EnumGroup<CoordinateModeEnum>::SmartPtr outputModeGroup = cl.AddEnum( "output-coordinates", &outputMode, "Coordinate specification mode for program output." );
    outputModeGroup->AddSwitch( Key( "absolute" ), COORDINATES_ABSOLUTE, "Use absolute volume coordinates. For each dimension, the valid range is [0,FOV]." );
    outputModeGroup->AddSwitch( Key( "indexed" ), COORDINATES_INDEXED, "Use grid indexes to specify coordinates. For each dimension, the valid value range is [0,Dims-1]." );
    outputModeGroup->AddSwitch( Key( "relative" ), COORDINATES_RELATIVE, "Use relative volume coordinates. For each dimension, the valid range is [0,1]." );
    outputModeGroup->AddSwitch( Key( "physical" ), COORDINATES_PHYSICAL, "Use physical volume coordinates. "
				"Each given location is transformed into image coordinates via the inverse of the images's index-to-physical space matrix." );

    cl.AddOption( Key( "radius" ), &radiusStr, "Radius of the search region in pixels (specified either as triple \"rX,rY,rZ\", or a single value, \"rXYZ\"). "
		  "The region searched is [2*rX+1,2*rY+1,2*rZ+1] pixels large, centered at the input location (but cropped at the image boundary)." );
    cl.AddSwitch( Key( "no-reorient" ), &readOrientation, static_cast<const char*>( NULL ), "Disable image reorientation into RAS alignment." );

    cl.AddParameter( &inputImagePath, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr volume;
  if ( readOrientation )
    volume = cmtk::VolumeIO::ReadOriented( inputImagePath, readOrientation );
  else
    volume = cmtk::VolumeIO::Read( inputImagePath );
      
  if ( ! volume || ! volume->GetData() ) 
    {
    cmtk::StdErr << "ERROR: could not read image " << inputImagePath << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::FixedVector<3,int> radius;
  if ( 3 != sscanf( radiusStr, "%6d,%6d,%6d", &radius[0], &radius[1], &radius[2] ) )
    {
    radius[0] = radius[1] = radius[2] = atof( radiusStr );
    }
  
  cmtk::UniformVolume::SpaceVectorType v;
  std::string restOfLine;

  while ( ! std::cin.eof() )
    {
    std::cin >> v;
    std::getline( std::cin, restOfLine );

    if ( std::cin.eof() )
      break;

    switch ( inputMode )
      {
      case COORDINATES_INDEXED:
	// nothing to do - lookup will be done by absolute coordinate
	break;
      case COORDINATES_ABSOLUTE:
	// index is absolute image coordinate divided (component-wise) by pixel size
	v /= volume->m_Delta;
	break;
      case COORDINATES_RELATIVE:
	// absolute image coordinate is relative times volume size
	(v *= volume->m_Size) /= volume->m_Delta;
	break;
      case COORDINATES_PHYSICAL:
	// absolute image coordinate is physical transformed by inverse image-to-physical matrix
	(v *= volume->GetImageToPhysicalMatrix().GetInverse()) /= volume->m_Delta;
	break;
      }
    
    // convert to discrete index by rounding (not simple truncation)
    const cmtk::UniformVolume::IndexType ijk = v.AddScalar( 0.5 );

    cmtk::DataGrid::IndexType maxIndex = ijk;
    cmtk::Types::DataItem maxValue = volume->GetDataAt( volume->GetOffsetFromIndex( ijk ) );

    cmtk::DataGrid::IndexType probe = ijk;
    for ( probe[2] = std::max( 0, ijk[2]-radius[2] ); probe[2] < std::min( volume->m_Dims[2], ijk[2]+radius[2]+1 ); ++probe[2] )
      for ( probe[1] = std::max( 0, ijk[1]-radius[1] ); probe[1] < std::min( volume->m_Dims[1], ijk[1]+radius[1]+1 ); ++probe[1] )
	for ( probe[0] = std::max( 0, ijk[0]-radius[0] ); probe[0] < std::min( volume->m_Dims[0], ijk[0]+radius[0]+1 ); ++probe[0] )
	  {
	  const cmtk::Types::DataItem value = volume->GetDataAt( volume->GetOffsetFromIndex( probe ) );
	  if ( value > maxValue )
	    {
	    maxValue = value;
	    maxIndex = probe;
	    }
	  }
    
    v = maxIndex;
    switch ( outputMode )
      {
      case COORDINATES_INDEXED:
	// nothing to do - already indexed
	break;
      case COORDINATES_ABSOLUTE:
	// absolute image coordinate is index times pixel size
	v *= volume->m_Delta;
	break;
      case COORDINATES_RELATIVE:
	// absolute image coordinate is relative times volume size
	(v *= volume->m_Delta) /= volume->m_Size;
	break;
      case COORDINATES_PHYSICAL:
	// absolute image coordinate is obtained using image-to-physical matrix
	(v *= volume->m_Delta) *= volume->GetImageToPhysicalMatrix();
	break;
      }
    
    std::cout << v << restOfLine << "\n";
    }

  return 0;
}

#include "cmtkSafeMain" 

