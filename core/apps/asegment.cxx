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

#include <cmtkCommandLine.h>
#include <cmtkConsole.h>

#include <list>

#include <cmtkVolumeIO.h>
#include <cmtkAtlasSegmentation.h>

bool Verbose = false;

const char* RawImageName = NULL;
const char* AtlasImageName = NULL;
const char* AtlasLabelName = NULL;
const char* OutFileName = "segmentation.hdr";

std::list<const char*> AtlasList;

bool
ParseCommandLine( const int argc, const char* argv[] )
{
  cmtk::CommandLine  cl( argc, argv );

  cl.AddSwitch( cmtk::CommandLine::Key( 'v', "verbose" ), &Verbose, true, "Verbose mode." );
  cl.AddOption( cmtk::CommandLine::Key( 'o', "output" ), &OutFileName, "Output file name [format:]path" );

  if ( ! cl.Parse() ) return false;

  // get raw image name
  RawImageName = cl.GetNextOptional();
  
  // get at least atlas images and label field.
  const char* next = cl.GetNextOptional();
  while ( next ) 
    {
    AtlasList.push_back( next );
    next = cl.GetNextOptional();
    }
  
  if ( ! RawImageName ) 
    {
    cmtk::StdErr << "ERROR: no images given. Use -h switch for help.\n";
    return false;
    }
  
  if ( AtlasList.size() % 2 ) 
    {
    cmtk::StdErr << "ERROR: must have alternating images and label fields for all atlases.\n";
    return false;
    }
  
  if ( AtlasList.size() < 2 ) 
    {
    cmtk::StdErr << "ERROR: must have at least one atlas image and labelmap.\n";
    return false;
    }
  
  return true;
}

int
main( const int argc, const char* argv[] )
{
  if ( ! ParseCommandLine( argc, argv ) ) 
    {
    exit( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr rawImage( cmtk::VolumeIO::ReadOriented( RawImageName, Verbose ) );
  if ( !rawImage ) 
    {
    cmtk::StdErr << "ERROR: could not read unsegmented image " << RawImageName << "\n";
    exit( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr atlasImg( cmtk::VolumeIO::ReadOriented( AtlasImageName, Verbose ) );
  if ( !atlasImg ) 
    {
    cmtk::StdErr << "ERROR: could not read atlas image " << AtlasImageName << "\n";
    exit( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr atlasLbl( cmtk::VolumeIO::ReadOriented( AtlasLabelName, Verbose ) );
  if ( !atlasLbl ) 
    {
    cmtk::StdErr << "ERROR: could not read atlas labels " << AtlasLabelName << "\n";
    exit( 1 );
    }
    
  if ( Verbose ) 
    {
    cmtk::StdErr << "Computing combined segmentation...\n";
    }

  cmtk::AtlasSegmentation segment( rawImage, atlasImg, atlasLbl, Verbose );
  cmtk::VolumeIO::Write( segment.GetLabelMap(), OutFileName );

  return 0;
}
