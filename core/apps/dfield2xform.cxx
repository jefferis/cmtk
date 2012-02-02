/*
//
//  Copyright 1997-2011 Torsten Rohlfing
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
#include <System/cmtkProgress.h>

#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkDeformationField.h>
#include <Base/cmtkFitSplineWarpToDeformationField.h>

#include <IO/cmtkXformIO.h>

const char* InputPath = NULL;
const char *OutputPath = NULL;

cmtk::Types::Coordinate GridSpacing = 0;

bool Absolute = false;

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
    cl.AddSwitch( Key( "absolute" ), &Absolute, true, "Input is absolute transformation field, x->u(x)." );
    cl.AddSwitch( Key( "relative" ), &Absolute, false, "Input is relative deformation field, x->x+u(x)." );
    cl.EndGroup();

    cl.BeginGroup( "Output", "Output Options" );
    cl.AddOption( Key( "grid-spacing" ), &GridSpacing, "Control point grid spacing of the output B-spline transformation." );
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

  cmtk::DeformationField::SmartPtr dfield = cmtk::DeformationField::SmartPtr::DynamicCastFrom( cmtk::XformIO::Read( InputPath ) );
  
  cmtk::FitSplineWarpToDeformationField fitSpline( dfield, Absolute );
  cmtk::SplineWarpXform::SmartPtr splineWarp = fitSpline.Fit( GridSpacing );

  cmtk::XformIO::Write( splineWarp, OutputPath );

  return 0;
}

#include "cmtkSafeMain"
