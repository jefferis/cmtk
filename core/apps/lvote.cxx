/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <IO/cmtkVolumeIO.h>

#include <vector>
#include <string>

#include <Unstable/cmtkLabelCombinationLocalVoting.h>

int
doMain
( const int argc, const char *argv[] )
{
  const char* targetImagePath = NULL;
  std::vector<std::string> atlasImagesLabels;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Local voting." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool combines multiple segmentations fro co-registered and reformatted atlases." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] targetImage atlasIntensity1 atlasLabels1 [atlasIntensity2 atlasLabels2 [...]]" );

    typedef cmtk::CommandLine::Key Key;

    cl.AddParameter( &targetImagePath, "TargetImage", "Target image path. This is the image to be segmented." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameterVector( &atlasImagesLabels, "AtlasImagesLabels", "List of reformatted atlas intensity and label images. This must be an even number of paths, where the first path within each pair is the intensity channel of"
			   "an atlas, and the second a label map channel of the same atlas, each reformatted into the space of the target image via an appropriate registration.");

    cl.Parse( argc, argv );

    if ( atlasImagesLabels.size() < 2 )
      {
      throw cmtk::CommandLine::Exception( "List of atlas intensity and label images must have at least two entries (one image and one label map file)" );
      }
    if ( atlasImagesLabels.size() % 2 )
      {
      throw cmtk::CommandLine::Exception( "List of atlas intensity and label images must have an even number of entries (one image and one label map file per atlas)" );
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartConstPtr targetImage = cmtk::VolumeIO::Read( targetImagePath );

  cmtk::LabelCombinationLocalVoting lvote( targetImage );
  for ( size_t atlasIdx = 0; atlasIdx < atlasImagesLabels.size(); )
    {
    lvote.AddAtlas( cmtk::VolumeIO::Read( atlasImagesLabels[atlasIdx++].c_str() ), 
		    cmtk::VolumeIO::Read( atlasImagesLabels[atlasIdx++].c_str() ) );
    }

  return 0;
}

#include "cmtkSafeMain"
