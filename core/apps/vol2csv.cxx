/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <IO/cmtkVolumeIO.h>

#include <vector>
#include <string>

int
doMain
( const int argc, const char *argv[] )
{
  const char* regionsImagePath = NULL;
  const char* labelsFilePath = NULL;
  const char* pxvolImagePath = NULL;
  const char* outputFilePath = NULL;

  std::vector<std::string> densityImagePaths;

  cmtk::Types::DataItem normalizeDensities = 1.0;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Compute regional volumes and write to CSV file." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool computes the volumes of regions in a label image. "
		       "It optionally accepts density maps (e.g., for different tissues) and computes and prints the per-region content for each. "
		       "Also, the tool can accept an optional 'pixel volume' map to account for local pixel volume variations, e.g., due to spatial distortion." );

    cl.AddParameter( &regionsImagePath, "RegionsImage", "Image of labeled regions." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameterVector( &densityImagePaths, "DensityImages", "List of density images. For each image given here, the total density per region is computed for each label in the regions image.")
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OPTIONAL );

    typedef cmtk::CommandLine::Key Key;

    cl.BeginGroup( "input", "Input Options" );
    cl.AddOption( Key( "normalize-densities" ), &normalizeDensities, "Optional normalization factor for density images. Typically, the values in the density images should be in the range 0..1, but often such images are scaled to "
		  "different ranges to accomodate storage as integers. If, for example, densities are stored as values 0..255, set this paramater to 255." );
    cl.AddOption( Key( "labels-file" ), &labelsFilePath, "If provided, this text file contains names for all labels in the regions image. These names are then used to label the rows of the CSV output." );
    cl.AddOption( Key( "pixel-volumes-image" ), &pxvolImagePath, "If provided, this volume contains scale factors for the volume of each pixel. This is typically the Jacobian determinant map of a spatial unwarping deformation." );
    cl.EndGroup();

    cl.BeginGroup( "output", "Output Options" );
    cl.AddOption( Key( 'o', "output" ), &outputFilePath, "If provided, program output is written to this file. If not provided, output is written to the STDOUT stream." );
    cl.EndGroup();

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartConstPtr regionsImage( cmtk::VolumeIO::Read( regionsImagePath ) );
  if ( ! regionsImage )
    {
    cmtk::StdErr << "ERROR: could not read regions image " << regionsImagePath << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr pxvolImage;
  if ( pxvolImagePath )
    {
    pxvolImage =  cmtk::VolumeIO::Read( pxvolImagePath );

    if ( ! pxvolImage )
      {
      cmtk::StdErr << "ERROR: could not read pixel volume image " << pxvolImagePath << "\n";
      throw cmtk::ExitException( 1 );
      }

    if ( ! regionsImage->GridMatches( *pxvolImage ) )
      {
      cmtk::StdErr << "ERROR: grid of pixel volume image " << pxvolImagePath << " does not match that of the regions image.\n";
      throw cmtk::ExitException( 1 );
      }
    }

  std::vector<cmtk::UniformVolume::SmartConstPtr> densityImages;
  for ( size_t idx = 0; idx < densityImagePaths.size(); ++idx )
    {
    cmtk::UniformVolume::SmartPtr nextImage( cmtk::VolumeIO::Read( densityImagePaths[idx].c_str() ) );
    if ( ! nextImage )
      {
      cmtk::StdErr << "ERROR: could not read density image " << densityImagePaths[idx] << "\n";
      throw cmtk::ExitException( 1 );
      }

    if ( ! regionsImage->GridMatches( *nextImage ) )
      {
      cmtk::StdErr << "ERROR: grid of density image " << densityImagePaths[idx] << " does not match that of the regions image.\n";
      throw cmtk::ExitException( 1 );
      }

    if ( normalizeDensities != 1.0 )
      nextImage->GetData()->Rescale( 1.0 / normalizeDensities );

    densityImages.push_back( nextImage );
    }

  const cmtk::Types::Coordinate pixelVolumeRegionsImage = regionsImage->m_Delta.Product();
  if ( pxvolImage )
    pxvolImage->GetData()->Rescale( pixelVolumeRegionsImage );

  const size_t maxLabel = std::max( 1, static_cast<int>( regionsImage->GetData()->GetRange().m_UpperBound ) );
  std::vector<cmtk::Types::Coordinate> regionVolumes( 1+maxLabel, 0.0 );

  std::vector< std::vector<cmtk::Types::Coordinate> > regionDensities( densityImages.size() );
  for ( size_t midx = 0; midx < densityImages.size(); ++midx )
    {
    regionDensities[midx].resize( 1+maxLabel, 0.0 );
    }
      
  cmtk::Types::Coordinate pixelVolume = 0.0;
  for ( size_t px = 0; px < regionsImage->GetNumberOfPixels(); ++px )
    {
    const size_t label = std::min<size_t>( maxLabel, std::max( 0, static_cast<int>( regionsImage->GetDataAt( px ) ) ) );

    // get size for this pixel depending on whether we have a per-pixel map or not.
    if ( pxvolImage )
      pixelVolume = pxvolImage->GetDataAt( px );
    else
      pixelVolume = pixelVolumeRegionsImage;

    regionVolumes[label] += pixelVolume;
    
    for ( size_t midx = 0; midx < densityImages.size(); ++midx )
      {
      regionDensities[midx][label] += pixelVolume * densityImages[midx]->GetDataAt( px );
      }
    }
  
  for ( size_t label = 0; label <= maxLabel; ++label )
    {
    std::cout << label << "," << regionVolumes[label];
    for ( size_t midx = 0; midx < densityImages.size(); ++midx )
      {
      std::cout << "," << regionDensities[midx][label];
      }
    std::cout << "\n";
    }

  return 0;
}

#include "cmtkSafeMain"
