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

#include <Base/cmtkWarpXform.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkFitAffineToWarpXform.h>

#include <IO/cmtkXformIO.h>

std::string InputPath;
std::string OutputPath;

int
doMain ( const int argc, const char *argv[] ) 
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Fit Affine Transformation to Nonrigid Transformation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Fit a linear affine transformation to a nonrigid transformation, either a B-spline free-form deformation or a non-parametric deformation field." );
    
    cl.AddParameter( &InputPath, "InputDField", "Input transformation." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );  
    cl.AddParameter( &OutputPath, "OutputXform", "Path for output fitted affine transformation." )->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  cmtk::WarpXform::SmartPtr warpXform = cmtk::WarpXform::SmartPtr::DynamicCastFrom( cmtk::XformIO::Read( InputPath ) );
  
  cmtk::FitAffineToWarpXform fitAffine( warpXform );
  cmtk::XformIO::Write( fitAffine.Fit(), OutputPath );
  
  return 0;
}

#include "cmtkSafeMain"
