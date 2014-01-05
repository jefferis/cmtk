/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>

#include <IO/cmtkXformIO.h>
#include <IO/cmtkVolumeIO.h>

#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkDeformationField.h>

std::string inXformPath;
std::string outXformPath;

float Fractional = -1;
bool DeformationOnly = false;

int
doMain ( const int argc, const char* argv[] ) 
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Convert nonrigd transformations." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool converts nonrigid B-spline free-format deformation coordinate transformations between different representations (e.g., absolute vs. relative vectors). Also creates fractional transformations." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( 'f', "fractional" ), &Fractional, "Write fractional deformation. Range: 0=affine to 1=full nonrigid; Default: 1" );
    cl.AddSwitch( Key( 'd', "deformation-only" ), &DeformationOnly, true, "Write only deformation part of transformation (minus global affine component)" );

    cl.AddParameter( &inXformPath, "InputPath", "Input transformation path" )->SetProperties( cmtk::CommandLine::PROPS_XFORM );
    cl.AddParameter( &outXformPath, "OutputPath", "Output transformation path" )->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( inXformPath ) );

  cmtk::SplineWarpXform::SmartPtr splineWarpXform = cmtk::SplineWarpXform::SmartPtr::DynamicCastFrom( xform );
  if ( splineWarpXform )
    {
    cmtk::AffineXform::SmartPtr initialAffine = splineWarpXform->GetInitialAffineXform();

    try
      {
      splineWarpXform->ReplaceInitialAffine();
      }
    catch ( cmtk::AffineXform::MatrixType::SingularMatrixException& ex )
      {
      cmtk::StdErr << "ERROR: singular matrix encountered in call to cmtk::WarpXform::ReplaceInitialAffine()\n";
      throw cmtk::ExitException( 1 );
      }
    
    if ( (Fractional >= 0) && (Fractional <= 1) )
      {
      const size_t numberOfControlPoints = splineWarpXform->GetNumberOfControlPoints();
      for ( size_t idx = 0; idx < numberOfControlPoints; ++idx )
	{
	cmtk::Vector3D v0 = splineWarpXform->GetOriginalControlPointPositionByOffset( idx );
	cmtk::Vector3D v1 = splineWarpXform->GetShiftedControlPointPositionByOffset( idx );

	((v1 -= v0) *= Fractional) += v0;

	splineWarpXform->SetShiftedControlPointPositionByOffset( v1, idx );
	}
      }
    
    if ( ! DeformationOnly )
      {
      try
	{
	splineWarpXform->ReplaceInitialAffine( initialAffine );
	}
      catch ( cmtk::AffineXform::MatrixType::SingularMatrixException& ex )
	{
	cmtk::StdErr << "ERROR: singular matrix encountered in call to cmtk::WarpXform::ReplaceInitialAffine()\n";
	throw cmtk::ExitException( 1 );
	}
      }
    }
  
  cmtk::XformIO::Write( xform, outXformPath );
  return 0;
}

#include "cmtkSafeMain"
