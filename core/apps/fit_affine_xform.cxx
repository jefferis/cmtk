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
//  $Revision: 4432 $
//
//  $LastChangedDate: 2012-06-13 13:33:51 -0700 (Wed, 13 Jun 2012) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkConsole.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkProgress.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkFitAffineToXformList.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkXformListIO.h>

int
doMain ( const int argc, const char *argv[] ) 
{
  std::string inputImagePath;
  std::string outputPath;
  
  bool fitRigid = false;

  cmtk::Types::Coordinate inversionTolerance = 1e-8;
  std::vector<std::string> inputXformPaths;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Fit Single Affine Transformation to Concatenated List" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Fit a linear affine transformation to a list of concatenated, optionally inverted, transformations." );
    
    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Input", "Input Options" );
    cl.AddOption( Key( "inversion-tolerance" ), &inversionTolerance, "Numerical tolerance of B-spline inversion in mm. Smaller values will lead to more accurate inversion, but may increase failure rate." );
    cl.EndGroup();

    cl.BeginGroup( "Fitting", "Fitting Options" );
    cl.AddSwitch( Key( "rigid" ), &fitRigid, true, "Fit rigid transformation (rotation and translation only) using SVD." );
    cl.AddSwitch( Key( "affine" ), &fitRigid, false, "Fit full affine transformation (rotation, translation, scales, shears) using pseudoinverse." );
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

  cmtk::XformList xformList = cmtk::XformListIO::MakeFromStringList( inputXformPaths );
  xformList.SetEpsilon( inversionTolerance );

  cmtk::UniformVolume::SmartPtr imageGrid( cmtk::VolumeIO::ReadGridOriented( inputImagePath ) );
  
  cmtk::FitAffineToXformList fit( *imageGrid, xformList );
  try
    {
    cmtk::AffineXform::SmartPtr xform = fit.Fit( fitRigid );
    cmtk::XformIO::Write( xform, outputPath );
    }
  catch ( cmtk::AffineXform::MatrixType::SingularMatrixException& ex )
    {
    cmtk::StdErr << "ERROR: singular matrix encountered in cmtk::FitAffineToXformList::Fit()\n";
    return 1;
    }

  return 0;
}

#include "cmtkSafeMain"
