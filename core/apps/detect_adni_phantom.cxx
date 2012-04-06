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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#include <cmtkconfig.h>

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>

#include <IO/cmtkVolumeIO.h>

#include <Segmentation/cmtkDetectPhantomMagphanEMR051.h>

#include <vector>

int
doMain( const int argc, const char* argv[] )
{
  const char* inputPath = NULL;
  const char* outputPath = NULL;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Detect ADNI phantom landmarks in phantom image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool detects the locations of all spherical landmarks in a 3D image of the Magphan EMR051 structural imaging phantom (a.k.a. ADNI Phantom)." );

    typedef cmtk::CommandLine::Key Key;
    
    cl.AddParameter( &inputPath, "InputImage", "Input image path. This is the image in which spheres are detected." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outputPath, "OutputImage", "Output image path. This image contains the magnitude filter response of the input image with respect to a matched, bipolar spherical correlation filter." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( inputPath ) );

  cmtk::DetectPhantomMagphanEMR051 detectionFilter( volume );
  std::vector<cmtk::UniformVolume::SpaceVectorType> landmarks = detectionFilter.GetLandmarks();

  for ( size_t i = 0; i < 7; ++i )
    {
    cmtk::StdOut << landmarks[i] << "\n";
    }

//  cmtk::VolumeIO::Write( *volume, outputPath );

  return 0;
}

#include "cmtkSafeMain"
