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

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>

#include <cmtkVolumeIO.h>
#include <cmtkDataGrid.h>
#include <cmtkMathUtil.h>
#include <cmtkUniformVolume.h>

#include <algorithm>

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

cmtk::Types::Coordinate levelsetThreshold = 1.0;

const char* inFile = NULL;
const char* outFile = NULL;

int
main( int argc, char* argv[] )
{
  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Levelset segmentation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Levelset-type segmentation of foreground/background using minimum regional variance energy" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] inputImage outputImage" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'n', "iterations" ), &numberOfIterations, "Maximum number of iterations [default: 100]" );
    cl.AddSwitch( Key( 'f', "force-iterations" ), &forceIterations, true, "Force given number of iterations, even when convergence has been detected" );

    cl.AddOption( Key( 's', "filter-sigma" ), &filterSigma, "Gaussian filter sigma [default: 2.0 mm]" );
    cl.AddOption( Key( 'd', "delta" ), &delta, "Time constant for levelset evolution; must be > 0; larger is faster [default: 0.1]" );
    cl.AddOption( Key( 't', "levelset-threshold" ), &levelsetThreshold, "Levelset threshold: levelset function is truncated at +/- this value [default: 1.0]" );

    cl.Parse();

    inFile = cl.GetNext();
    outFile = cl.GetNext();
    }
  catch ( cmtk::CommandLine::Exception ex )
    {
    cmtk::StdErr << ex;
    exit( 1 );
    }

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( inFile, verbose ) );
  const size_t numberOfPixels = volume->GetNumberOfPixels();

  cmtk::UniformVolume::SmartPtr levelset( volume->CloneGrid() );
  levelset->CreateDataArray( cmtk::TYPE_FLOAT );
  levelset->GetData()->Fill( -1.0 );
  cmtk::Vector3D center( volume->GetDims(0)/2, volume->GetDims(1)/2, volume->GetDims(2)/2 );
  levelset->DrawSphere( center, (levelset->GetDims(0)+levelset->GetDims(1)+levelset->GetDims(2))/6, 1.0 );

  size_t nInsideOld = 0, nInside = 1;

  for ( int it = 0; (it < numberOfIterations) && ((nInside!=nInsideOld) || forceIterations); ++it )
    {
    nInsideOld = nInside;
    nInside = 0;
    double insideSum = 0, outsideSum = 0;

    levelset->ApplyGaussFilter( filterSigma );
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
    
    const float mInside = insideSum / nInside;
    const float mOutside = outsideSum / nOutside;

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
  
  cmtk::VolumeIO::Write( levelset, outFile, verbose );

  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace levelset
} // namespace apps
} // namespace cmtk
#endif

