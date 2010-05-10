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

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>
#include <cmtkProgressConsole.h>

#include <cmtkVolumeIO.h>
#include <cmtkMathUtil.h>
#include <cmtkUniformVolume.h>
#include <cmtkUniformVolumePainter.h>
#include <cmtkUniformVolumeFilter.h>

#include <algorithm>

#ifdef CMTK_USE_SQLITE
#  include <cmtkImageXformDB.h>
#endif

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
  const size_t numberOfPixels = volume->GetNumberOfPixels();

  cmtk::UniformVolume::SmartPtr levelset( volume->CloneGrid() );
  levelset->CreateDataArray( cmtk::TYPE_FLOAT );
  levelset->GetData()->Fill( -1.0 );

  cmtk::Vector3D center( volume->GetDims()[0]/2, volume->GetDims()[1]/2, volume->GetDims()[2]/2 );

  cmtk::UniformVolumePainter painter( levelset );
  painter.DrawSphere( center, (levelset->GetDims()[0]+levelset->GetDims()[1]+levelset->GetDims()[2])/6, 1.0 );

  size_t nInsideOld = 0, nInside = 1;

  cmtk::Progress::Begin( 0, numberOfIterations, 1, "Levelset Evolution" );
  for ( int it = 0; (it < numberOfIterations) && ((nInside!=nInsideOld) || forceIterations); ++it )
    {
    cmtk::Progress::SetProgress( it );

    nInsideOld = nInside;
    nInside = 0;
    cmtk::Types::DataItem insideSum = 0, outsideSum = 0;

    levelset->SetData( cmtk::UniformVolumeFilter( levelset ).GetDataGaussFiltered( filterSigma ) );
#pragma omp parallel for reduction(+:nInside) reduction(+:insideSum) reduction(+:outsideSum)
    for ( size_t n = 0; n < numberOfPixels; ++n )
      {
      if ( levelset->GetDataAt( n ) > 0 )
	{
	insideSum += volume->GetDataAt( n );
	++nInside;
	}
      else
	outsideSum += volume->GetDataAt( n );
      }

    const size_t nOutside = numberOfPixels - nInside;
    const cmtk::Types::DataItem ratioInOut = 1.0 * nInside / nOutside;
    
    const cmtk::Types::DataItem mInside = insideSum / nInside;
    const cmtk::Types::DataItem mOutside = outsideSum / nOutside;

    if ( verbose )
      std::cerr << it << " IN: " << nInside << "  " << mInside << "  OUT: " << nOutside << "  " << mOutside << "\r";
    
#pragma omp parallel for
    for ( size_t n = 0; n < numberOfPixels; ++n )
      {
      const cmtk::Types::DataItem data = volume->GetDataAt( n );
      const cmtk::Types::DataItem zInside = fabs( mInside - data );
      const cmtk::Types::DataItem zOutside = fabs( mOutside - data );
      cmtk::Types::DataItem newLevel = levelset->GetDataAt( n );
      if ( zInside>zOutside )
	{
	newLevel -= delta * ratioInOut;
	}
      else
	{
	newLevel += delta / ratioInOut;
	}
      levelset->SetDataAt( std::min<cmtk::Types::DataItem>( levelsetThreshold, std::max<cmtk::Types::DataItem>( -levelsetThreshold, newLevel ) ), n );
      }
    }

  cmtk::Progress::Done();

  if ( binarize )
    {
    levelset->GetData()->Binarize( 0.5 );
    levelset->SetData( cmtk::TypedArray::SmartPtr( levelset->GetData()->Convert( cmtk::TYPE_BYTE ) ) );
    }
  
  cmtk::VolumeIO::Write( *levelset, outFile, verbose );

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

