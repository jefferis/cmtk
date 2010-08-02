/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "System/cmtkConsole.h"
#include "System/cmtkCommandLine.h"
#include "System/cmtkProgressConsole.h"

#include "IO/cmtkVolumeIO.h"

#include <algorithm>

#ifdef CMTK_USE_SQLITE
#  include "Registration/cmtkImageXformDB.h"
#endif

#include "Segmentation/cmtkSimpleLevelset.h"

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace levelset
{
#endif
bool verbose = false;

cmtk::Types::Coordinate filterSigma = 2.0;
cmtk::Types::Coordinate delta = 0.1;

int numberOfIterations = 100;
bool forceIterations = false;
bool binarize = false;

cmtk::Types::Coordinate levelsetThreshold = 1.0;

const char* inFile = NULL;
const char* outFile = NULL;

#ifdef CMTK_USE_SQLITE
const char* updateDB = NULL;
#endif

int
main( int argc, char* argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv, cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Levelset segmentation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Levelset-type segmentation of foreground/background using minimum regional variance energy" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Segmentation" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose mode" )->SetProperties( cmtk::CommandLine::PROPS_NOXML );

    cl.AddSwitch( Key( 'b', "binarize" ), &binarize, true, "Binarize levelset and write as byte mask, rather than write floating-point levelset function itself." );

    cl.BeginGroup( "Levelset Evolution Parameters", "These parameters of control the evolution of the levelset function" )->SetProperties( cmtk::CommandLine::PROPS_ADVANCED );
    cl.AddOption( Key( 'n', "iterations" ), &numberOfIterations, "Maximum number of iterations" );
    cl.AddSwitch( Key( 'f', "force-iterations" ), &forceIterations, true, "Force given number of iterations, even when convergence has been detected" );

    cl.AddOption( Key( 's', "filter-sigma" ), &filterSigma, "Gaussian filter sigma in world coordinate units (e.g., mm)" );
    cl.AddOption( Key( 'd', "delta" ), &delta, "Time constant for levelset evolution; must be > 0; larger is faster" );
    cl.AddOption( Key( 't', "levelset-threshold" ), &levelsetThreshold, "Levelset threshold: levelset function is truncated at +/- this value" );
    cl.EndGroup();

#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &updateDB, "Path to image/transformation database that should be updated with the newly created image." );
    cl.EndGroup();
#endif

    cl.AddParameter( &inFile, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outFile, "OutputImage", "Output image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_LABELS | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse();
    }
  catch ( const cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex;
    exit( 1 );
    }

  // Instantiate programm progress indicator.
  cmtk::ProgressConsole progressIndicator( "LevelsetSegmentation" );

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( inFile, verbose ) );

  cmtk::SimpleLevelset levelset( volume );
  levelset.SetFilterSigma( cmtk::Units::GaussianSigma( filterSigma ) );
  levelset.SetTimeDelta( delta );
  levelset.SetLevelsetThreshold( levelsetThreshold );

  levelset.InitializeCenteredSphere();
  levelset.Evolve( numberOfIterations, forceIterations );

  cmtk::VolumeIO::Write( *levelset.GetLevelset( binarize ), outFile, verbose );

#ifdef CMTK_USE_SQLITE
  if ( updateDB )
    {
    cmtk::ImageXformDB db( updateDB );
    db.AddImage( outFile, inFile );
    }
#endif

  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace levelset
} // namespace apps
} // namespace cmtk
#endif

