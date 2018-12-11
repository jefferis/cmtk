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

#include <IO/cmtkVolumeIO.h>

#include <Segmentation/cmtkSphereDetectionBipolarMatchedFilterFFT.h>
#include <Segmentation/cmtkSphereDetectionNormalizedBipolarMatchedFilterFFT.h>

int
doMain( const int argc, const char* argv[] )
{
  std::string inputPath;
  std::string outputPath;

  cmtk::Types::Coordinate sphereRadius = 1;
  int filterMargin = 2;
  bool normalized = false;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Detect spheres" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool detects spherical objects in three-dimensional images." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( 'r', "radius" ), &sphereRadius, "Radius of spheres to detect in physical length units (typically mm)." );
    cl.AddOption( Key( "filter-margin" ), &filterMargin, "Half of filter margin width in pixels. This determines the extent of the region around the surface or the sphere where the detection filter is non-zero." );
    cl.AddSwitch( Key( "normalized" ), &normalized, true, "Use intensity-normalized filter, effectively computing Normalized Cross Correlation between filter and image." );
    
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

  if ( normalized )
    {
    cmtk::SphereDetectionNormalizedBipolarMatchedFilterFFT detectionFilter( *volume );
    volume->SetData( detectionFilter.GetFilteredImageData( sphereRadius, filterMargin ) );
    }
  else
    {
    cmtk::SphereDetectionBipolarMatchedFilterFFT detectionFilter( *volume );
    volume->SetData( detectionFilter.GetFilteredImageData( sphereRadius, filterMargin ) );
    }

  cmtk::VolumeIO::Write( *volume, outputPath );

  return 0;
}

#include "cmtkSafeMain"
