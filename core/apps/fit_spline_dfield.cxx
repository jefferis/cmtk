/*
//
//  Copyright 1997-2011 Torsten Rohlfing
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

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkProgress.h>

#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkDeformationField.h>
#include <Base/cmtkFitAffineToWarpXform.h>
#include <Base/cmtkFitSplineWarpToDeformationField.h>

#include <IO/cmtkXformIO.h>

const char* InputPath = NULL;
const char *OutputPath = NULL;

const char* GridDims = NULL;
cmtk::Types::Coordinate GridSpacing = 0;
int Levels = 1;

bool AffineFirst = false;

bool Absolute = true;

int
doMain ( const int argc, const char *argv[] ) 
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Deformation Field to Transformation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Fit a parametric nonrigid transformation (B-spline free-form deformation) to a deformation field" );
    
    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Input", "Input Options" );
    cl.AddSwitch( Key( "absolute" ), &Absolute, true, "Input is an absolute transformation field [x maps to input(x)]" );
    cl.AddSwitch( Key( "relative" ), &Absolute, false, "Input is relative deformation field, e.g., a gradient or force field [x maps to x+input(x)]" );
    cl.EndGroup();

    cl.BeginGroup( "Fitting", "Fitting Options" );
    cl.AddOption( Key( "levels" ), &Levels, "Number of levels in the multi-level B-spline approximation procedure." );
    cl.AddSwitch( Key( "fit-affine-first" ), &AffineFirst, true, "Fit affine transformation first, then initialize spline with it." );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( "final-cp-spacing" ), &GridSpacing, "Final control point grid spacing of the output B-spline transformation." );
    cl.AddOption( Key( "final-cp-dims" ), &GridDims, "Final control point grid dimensions (i.e., number of controlpoints) of the output B-spline transformation. To be provided as 'dimX,dimY,dimZ'." );
    cl.EndGroup();

    cl.AddParameter( &InputPath, "InputDField", "Input deformation field." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );  
    cl.AddParameter( &OutputPath, "OutputXform", "Path for the output transformation." )->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
 
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  if ( GridDims && GridSpacing )
    {
    cmtk::StdErr << "ERROR: must specify either output spline control point spacing or grid dimensions, but not both.\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::DeformationField::SmartPtr dfield = cmtk::DeformationField::SmartPtr::DynamicCastFrom( cmtk::XformIO::Read( InputPath ) );
  
  cmtk::FitSplineWarpToDeformationField fitSpline( dfield, Absolute );
  cmtk::SplineWarpXform::SmartPtr splineWarp;

  cmtk::AffineXform::SmartPtr affineXform;
  if ( AffineFirst )
    {
    affineXform = cmtk::FitAffineToWarpXform( dfield ).Fit();
    }

  if ( GridSpacing )
    {
    splineWarp = fitSpline.Fit( GridSpacing, Levels, affineXform );
    }
  else
    {
    if ( GridDims )
      {
      double dims[3];
      if ( 3 != sscanf( GridDims, "%lf,%lf,%lf", &(dims[0]), &(dims[1]), &(dims[2]) ) )
	{
	cmtk::StdErr << "ERROR: grid dimensions must be specified as dimsX,dimsY,dimsZ\n";
	throw cmtk::ExitException( 1 );
	}
      
      splineWarp = fitSpline.Fit( cmtk::FixedVector<3,double>::FromPointer( dims ), Levels, affineXform );
      }
    else
      {
      cmtk::StdErr << "ERROR: must specify either output spline control point spacing or grid dimensions.\n";
      throw cmtk::ExitException( 1 );
      }
    }
  
  cmtk::XformIO::Write( splineWarp, OutputPath );

  return 0;
}

#include "cmtkSafeMain"
