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
#include <cmtkProgressConsole.h>

#include <list>

#include <cmtkVolumeIO.h>
#include <cmtkAtlasSegmentation.h>

int
main( const int argc, const char* argv[] )
{
  bool verbose = false;
  bool fast = false;
  
  const char* targetImageName = NULL;
  const char* atlasImageName = NULL;
  const char* atlasLabelName = NULL;
  const char* outImageName = NULL;

  try
    {
    cmtk::CommandLine  cl( argc, argv, cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Atlas-based segmentation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Register a target image to an atlas, using affine followed by nonrigid B-spline registration, then reformat the atlas label map to the target image." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Segmentation" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose mode." );
    cl.AddSwitch( Key( 'f', "fast" ), &fast, true, "Fast mode." );
    
    cl.AddParameter( &targetImageName, "TargetImage", "Target image path. This is the image to be segmented." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &atlasImageName, "AtlasImage", "Atlas image path. This is the structural channel (e.g., T1-weighted MRI) of the atlas to be used." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &atlasLabelName, "AtlasLabels", "Atlas label image path. This is the label map to be reformatted to the target image." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_LABELS );
    cl.AddParameter( &outImageName, "OutputImage", "Output image path. This is where the reformatted label map is written." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_LABELS | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse();
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }
  
  // Instantiate programm progress indicator.
  cmtk::ProgressConsole progressIndicator( "Atlas-based Segmentation" );

  cmtk::UniformVolume::SmartPtr targetImg( cmtk::VolumeIO::ReadOriented( targetImageName, verbose ) );
  if ( !targetImg ) 
    {
    cmtk::StdErr << "ERROR: could not read target image " << targetImageName << "\n";
    exit( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr atlasImg( cmtk::VolumeIO::ReadOriented( atlasImageName, verbose ) );
  if ( !atlasImg ) 
    {
    cmtk::StdErr << "ERROR: could not read atlas image " << atlasImageName << "\n";
    exit( 1 );
    }
  
  cmtk::UniformVolume::SmartPtr atlasLbl( cmtk::VolumeIO::ReadOriented( atlasLabelName, verbose ) );
  if ( !atlasLbl ) 
    {
    cmtk::StdErr << "ERROR: could not read atlas labels " << atlasLabelName << "\n";
    exit( 1 );
    }
    
  cmtk::AtlasSegmentation segment( targetImg, atlasImg, atlasLbl );
  segment.SetVerbose( verbose );
  segment.SetFast( fast );

  cmtk::VolumeIO::Write( segment.GetLabelMap(), outImageName );

  return 0;
}
