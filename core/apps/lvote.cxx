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
//  $Revision: 2765 $
//
//  $LastChangedDate: 2011-01-18 16:17:56 -0800 (Tue, 18 Jan 2011) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>
#include <System/cmtkConsole.h>

#include <vector>
#include <string>

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
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] targetImage atlasImage1 atlasLabels1 [atlasImage2 atlasLabels2 [...]]" );

    typedef cmtk::CommandLine::Key Key;

    cl.AddParameter( &targetImagePath, "TargetImage", "Target image path. This is the image to be segmented." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameterVector( &atlasImagesLabels, "AtlasImagesLabels", "List of reformatted atlas images and label maps. This must be an even number of paths, where the first path within each pair is the image channel of"
			   "an atlas, and the second a label map channel of the same atlas, each reformatted into the space of the target image via an appropriate registration.");

    cl.Parse( argc, argv );

    if ( atlasImagesLabels.size() < 2 )
      {
      throw cmtk::CommandLine::Exception( "List of atlas images and label maps must have at least two entries (one image and one label map file)" );
      }
    if ( atlasImagesLabels.size() % 2 )
      {
      throw cmtk::CommandLine::Exception( "List of atlas images and label maps must have an even number of entries (one image and one label map file per atlas)" );
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  return 0;
}

#include "cmtkSafeMain"
