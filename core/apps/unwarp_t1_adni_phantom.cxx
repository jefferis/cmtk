/*
//
//  Copyright 1997-2010 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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
#include <System/cmtkDebugOutput.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>

#include <Segmentation/cmtkDetectPhantomMagphanEMR051.h>

#include <Base/cmtkLandmarkList.h>
#include <Base/cmtkLandmarkPairList.h>
#include <Base/cmtkFitAffineToLandmarks.h>
#include <Base/cmtkFitSplineWarpToLandmarks.h>

#include <vector>

int
doMain( const int argc, const char* argv[] )
{
  const char* inputPhantomPath = NULL;
  const char* inputImagePath = NULL;

  const char* gridDims = NULL;
  cmtk::Types::Coordinate gridSpacing = 0;
  int levels = 1;  
  bool affineFirst = false;

  const char* outputXform = NULL;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Unwarp T1-weighted MR image using ADNI phantom image" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool detects the locations of all spherical landmarks in a T1-weighted MR image of the Magphan EMR051 structural imaging phantom (a.k.a. ADNI Phantom) and computes a "
		       "B-spline free-form deformation to unwarp another T1-weighted image acquired on the same scanner." );

    typedef cmtk::CommandLine::Key Key;    
    cl.BeginGroup( "Fitting", "Fitting Options" );
    cl.AddOption( Key( "levels" ), &levels, "Number of levels in the multi-level B-spline approximation procedure." );
    cl.AddSwitch( Key( "fit-affine-first" ), &affineFirst, true, "Fit affine transformation first, then initialize spline with it." );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( "final-cp-spacing" ), &gridSpacing, "Final control point grid spacing of the output B-spline transformation." );
    cl.AddOption( Key( "final-cp-dims" ), &gridDims, "Final control point grid dimensions (i.e., number of controlpoints) of the output B-spline transformation. To be provided as 'dimX,dimY,dimZ'." );
    cl.EndGroup();

    cl.AddParameter( &inputPhantomPath, "InputPhantom", "Phantom image path. This is the image of the ADNI phantom, which is used to compute the unwarping transformation." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &inputImagePath, "InputPhantom", "Input image path. This is the image that is unwarped. It is important that this image be acquired on the same scanner (not only the same model but the very machine) "
		     "on which the phantom image was also acquired, preferably in close temporal proximity. Also, both this and the phantom image must share and specify the same physical image coordinates, i.e., only images in "
		     "NIFTI or NRRD format can be used." )
      ->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outputXform, "OutputXform", "Output transformation path. This is the affine phantom-to-image coordinate transformation fitted to the detected landmark spheres." )
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  if ( gridDims && gridSpacing )
    {
    cmtk::StdErr << "ERROR: must specify either output spline control point spacing or grid dimensions, but not both.\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartConstPtr phantomImage( cmtk::VolumeIO::ReadOriented( inputPhantomPath ) );
  cmtk::UniformVolume::SmartConstPtr unwarpImage( cmtk::VolumeIO::ReadOriented( inputImagePath ) );

  cmtk::DetectPhantomMagphanEMR051 detectionFilter( phantomImage );
  cmtk::LandmarkList expectedLandmarks = detectionFilter.GetExpectedLandmarks();
  cmtk::LandmarkList actualLandmarks = detectionFilter.GetDetectedLandmarks();

  cmtk::LandmarkPairList pairList( expectedLandmarks, actualLandmarks );
  cmtk::DebugOutput( 2 ) << "INFO: detected and matched " << pairList.size() << " out of " << expectedLandmarks.size() << " expected landmarks.\n";

  cmtk::SplineWarpXform::SmartConstPtr splineWarp;
  cmtk::AffineXform::SmartConstPtr affineXform;
  if ( affineFirst )
    {
    affineXform = cmtk::FitAffineToLandmarks( pairList ).GetAffineXform();
    }
  
  if ( gridSpacing )
    {
    splineWarp = cmtk::FitSplineWarpToLandmarks( pairList ).Fit( unwarpImage->Size, gridSpacing, levels, affineXform.GetPtr() );
    }
  else
    {
    if ( gridDims )
      {
      double dims[3];
      if ( 3 != sscanf( gridDims, "%lf,%lf,%lf", &(dims[0]), &(dims[1]), &(dims[2]) ) )
	{
	cmtk::StdErr << "ERROR: grid dimensions must be specified as dimsX,dimsY,dimsZ\n";
	throw cmtk::ExitException( 1 );
	}
      
      splineWarp = cmtk::FitSplineWarpToLandmarks( pairList ).Fit( unwarpImage->Size, cmtk::FixedVector<3,double>( dims ), levels, affineXform.GetPtr() );
      }
    else
      {
      cmtk::StdErr << "ERROR: must specify either output spline control point spacing or grid dimensions.\n";
      throw cmtk::ExitException( 1 );
      }
    }
  
  cmtk::XformIO::Write( splineWarp, outputXform );

  return 0;
}

#include "cmtkSafeMain"
