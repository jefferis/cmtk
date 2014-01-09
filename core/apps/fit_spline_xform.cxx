/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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
#include <System/cmtkProgress.h>

#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkFitSplineWarpToXformList.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkXformListIO.h>

int
doMain ( const int argc, const char *argv[] ) 
{
  std::string inputImagePath;
  std::string outputPath;
  
  std::string gridDims;
  cmtk::Types::Coordinate gridSpacing = 0;
  int levels = 1;
  
  cmtk::Types::Coordinate inversionTolerance = 1e-8;
  std::vector<std::string> inputXformPaths;

  bool affineFirst = false;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Fit Single B-Spline Transformation to Concatenated Transformation List" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Fit a parametric nonrigid transformation (B-spline free-form deformation) to a list of concatenated, optionally inverted, transformations." );
    
    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Fitting", "Fitting Options" );
    cl.AddOption( Key( "levels" ), &levels, "Number of levels in the multi-level B-spline approximation procedure." );
    cl.AddSwitch( Key( "fit-affine-first" ), &affineFirst, true, "Fit affine transformation first, then initialize spline with it." );
    cl.AddOption( Key( "final-cp-spacing" ), &gridSpacing, "Final control point grid spacing of the output B-spline transformation." );
    cl.AddOption( Key( "final-cp-dims" ), &gridDims, "Final control point grid dimensions (i.e., number of controlpoints) of the output B-spline transformation. To be provided as 'dimX,dimY,dimZ'." );
    cl.EndGroup();

    cl.BeginGroup( "Input", "Input Options" );
    cl.AddOption( Key( "inversion-tolerance" ), &inversionTolerance, "Numerical tolerance of B-spline inversion in mm. Smaller values will lead to more accurate inversion, but may increase failure rate." );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( 'o', "output" ), &outputPath, "Path for the output transformation." )->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    cl.EndGroup();

    cl.AddParameter( &inputImagePath, "InputImage", "Input image path. This image determines the discrete sampling grid where the target transformation is estimated and fitted." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameterVector( &inputXformPaths, "XformList", "List of concatenated transformations. Insert '--inverse' to use the inverse of the transformation listed next. "
			   "(If the first transformation in the sequence is inverted, then '--inverse' must be preceded by '--', i.e., use '-- --inverse xform.path')." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );
 
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  if ( !gridDims.empty() && gridSpacing )
    {
    cmtk::StdErr << "ERROR: must specify either output spline control point spacing or grid dimensions, but not both.\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::XformList xformList = cmtk::XformListIO::MakeFromStringList( inputXformPaths );
  xformList.SetEpsilon( inversionTolerance );

  cmtk::UniformVolume::SmartPtr imageGrid( cmtk::VolumeIO::ReadGridOriented( inputImagePath ) );
  
  cmtk::FitSplineWarpToXformList fitSpline( *imageGrid, xformList );
  cmtk::SplineWarpXform::SmartPtr splineWarp;
  
  if ( gridSpacing )
    {
    try
      {
      splineWarp = fitSpline.Fit( gridSpacing, levels, affineFirst );
      }
    catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
      {
      cmtk::StdErr << "ERROR: singular matrix encountered in cmtk::FitSplineWarpToXformList::Fit()\n";
      throw cmtk::ExitException( 1 );
      }   
    }
  else
    {
    if ( !gridDims.empty() )
      {
      double dims[3];
      if ( 3 != sscanf( gridDims.c_str(), "%20lf,%20lf,%20lf", &(dims[0]), &(dims[1]), &(dims[2]) ) )
	{
	cmtk::StdErr << "ERROR: grid dimensions must be specified as dimsX,dimsY,dimsZ\n";
	throw cmtk::ExitException( 1 );
	}

      try
	{
	splineWarp = fitSpline.Fit( cmtk::FixedVector<3,double>::FromPointer( dims ), levels, affineFirst );
	}
      catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& )
	{
	cmtk::StdErr << "ERROR: singular matrix encountered in cmtk::FitSplineWarpToXformList::Fit()\n";
	throw cmtk::ExitException( 1 );
	}
      }
    else
      {
      cmtk::StdErr << "ERROR: must specify either output spline control point spacing or grid dimensions.\n";
      throw cmtk::ExitException( 1 );
      }
    }
  
  cmtk::XformIO::Write( splineWarp, outputPath );

  return 0;
}

#include "cmtkSafeMain"
