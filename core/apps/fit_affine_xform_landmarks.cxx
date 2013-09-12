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
#include <System/cmtkDebugOutput.h>
#include <System/cmtkCommandLine.h>
#include <System/cmtkProgress.h>
#include <System/cmtkExitException.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkLandmarkPairList.h>
#include <Base/cmtkFitAffineToLandmarks.h>
#include <Base/cmtkFitRigidToLandmarks.h>

#include <IO/cmtkLandmarkListIO.h>
#include <IO/cmtkXformIO.h>

#include <fstream>

int
doMain ( const int argc, const char *argv[] ) 
{
  std::string lmSourcePath;
  std::string lmTargetPath;
  std::string outputPath;

  bool rigid = false;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Fit Affine Transformation to Landmarks" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Fit a linear affine transformation to a set of matched landmarks." );
    
    typedef cmtk::CommandLine::Key Key;
    cl.AddParameter( &lmSourcePath, "SourcePath", "Path to file with source-space landmarks." );
    cl.AddParameter( &lmTargetPath, "TargetPath", "Path to file with target-space landmarks." );

    cl.AddSwitch( Key( "rigid" ), &rigid, true, "Fit only a rigid, 6 degree-of-freedom transformation rather than full affine. This requires fewer landmarks and makes less trict assumptions about their location "
		  "(i.e., three non-collinear landmarks are sufficient to obtain a rigid transformation." );
    
    cl.AddParameter( &outputPath, "OutputXform", "Path for output fitted affine transformation." )->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT );
    
    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  std::ifstream lmSourceStream( lmSourcePath.c_str() );
  if ( !lmSourceStream.good() )
    {
    cmtk::StdErr << "ERROR: could not open source landmark file " << lmSourcePath << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::LandmarkList sourceLandmarks;
  lmSourceStream >> sourceLandmarks;

  cmtk::DebugOutput( 5 ) << "Read " << sourceLandmarks.size() << " source landmarks\n";

  std::ifstream lmTargetStream( lmTargetPath.c_str() );
  if ( !lmTargetStream.good() )
    {
    cmtk::StdErr << "ERROR: could not open target landmark file " << lmTargetPath << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::LandmarkList targetLandmarks;
  lmTargetStream >> targetLandmarks;

  cmtk::DebugOutput( 5 ) << "Read " << targetLandmarks.size() << " target landmarks\n";

  cmtk::LandmarkPairList landmarkPairs( sourceLandmarks, targetLandmarks );
  cmtk::DebugOutput( 5 ) << "Matched " << landmarkPairs.size() << " landmark pairs\n";

  cmtk::AffineXform::SmartConstPtr xform = rigid ? cmtk::FitRigidToLandmarks( landmarkPairs ).GetRigidXform() : cmtk::FitAffineToLandmarks( landmarkPairs ).GetAffineXform();
  cmtk::XformIO::Write( xform, outputPath );

  if ( cmtk::DebugOutput::GetGlobalLevel() >= 5 )
    {
    cmtk::DebugOutput( 5 ) << "INFO: RMS fitting residual = " << sqrt( xform->GetLandmarksMSD( landmarkPairs ) ) << "\n";
    }
  
  return 0;
}

#include "cmtkSafeMain"
