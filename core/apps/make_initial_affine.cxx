/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013 SRI International
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

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkMatrix4x4.h>
#include <Base/cmtkTransformChangeToSpaceAffine.h>

#include <Registration/cmtkMakeInitialAffineTransformation.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformIO.h>

#ifdef CMTK_USE_SQLITE
#  include <Registration/cmtkImageXformDB.h>
#endif

#include <stdlib.h>

int
doMain( const int argc, const char* argv[] )
{
  std::string referenceImagePath;
  std::string floatingImagePath;
  std::string outputXformPath;
  
  bool centerXform = false;

  bool writeXformNativeSpaces = false;
  if ( getenv( "Slicer3_HOME" ) != NULL )
    {
    writeXformNativeSpaces = true;
    }

  int mode = 0;

#ifdef CMTK_USE_SQLITE
  std::string updateDB;
#endif

  try
    {
    cmtk::CommandLine cl( cmtk::CommandLine::PROPS_XML );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Initialize affine transformation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Compute initial affine transformation by aligning centers of mass or principal axes" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_CATEG, "CMTK.Registration" );

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Transformation", "Transformation construction control" );
    cmtk::CommandLine::EnumGroup<int>::SmartPtr modeGroup = cl.AddEnum( "mode", &mode, "Mode selection for initialization" );
    modeGroup->AddSwitch( Key( "direction-vectors" ), 0, "Alignment based on image direction vectors" );
    modeGroup->AddSwitch( Key( "centers-of-mass" ), 1, "Alignment based on centers of mass (translation only)" );
    modeGroup->AddSwitch( Key( "principal-axes" ), 2, "Alignment based on principal axes" );
    modeGroup->AddSwitch( Key( "identity" ), -1, "Create only an identity transformation" );
    
    cl.AddSwitch( Key( 'C', "center-xform" ), &centerXform, true, "Set transformation center (for rotation, scale) to center of reference image." );
    cl.AddSwitch( Key( "native-space" ), &writeXformNativeSpaces, true, "Write transformation between the native image spaces, rather than in CMTK standard RAS space. This is the default when running this tool as a Slicer plugin." )
      ->SetProperties( cmtk::CommandLine::PROPS_NOXML );;
    cl.EndGroup();
    
#ifdef CMTK_USE_SQLITE
    cl.BeginGroup( "Database", "Image/Transformation Database" );
    cl.AddOption( Key( "db" ), &updateDB, "Path to image/transformation database that should be updated with the newly transformation from reference to floating image." );
    cl.EndGroup();
#endif

    cl.AddParameter( &referenceImagePath, "ReferenceImage", "Reference (fixed) image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &floatingImagePath, "FloatingImage", "Floating (moving) image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &outputXformPath, "OutputXform", "Output transformation path" )
      ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OUTPUT )
      ->SetAttribute( "reference", "FloatingImage" );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& ex )
    {
    cmtk::StdErr << ex << "\n";
    return 1;
    }

  cmtk::UniformVolume::SmartPtr referenceImage( cmtk::VolumeIO::ReadOriented( referenceImagePath ) );
  if ( ! referenceImage )
    {
    cmtk::StdErr << "ERROR: could not read image " << referenceImagePath << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::UniformVolume::SmartPtr floatingImage( cmtk::VolumeIO::ReadOriented( floatingImagePath ) );
  if ( ! floatingImage )
    {
    cmtk::StdErr << "ERROR: could not read image " << floatingImagePath << "\n";
    throw cmtk::ExitException( 1 );
    }
  
  cmtk::AffineXform::SmartPtr xform;
  switch ( mode )
    {
    case -1:
      xform = cmtk::AffineXform::SmartPtr( new cmtk::AffineXform );
      break;
    case 0:
      xform = cmtk::AffineXform::SmartPtr( cmtk::MakeInitialAffineTransformation::AlignDirectionVectors( *referenceImage, *floatingImage, centerXform ) );
      break;
    case 1:
      xform = cmtk::AffineXform::SmartPtr( cmtk::MakeInitialAffineTransformation::AlignCentersOfMass( *referenceImage, *floatingImage ) );
      break;
    case 2:
      xform = cmtk::AffineXform::SmartPtr( cmtk::MakeInitialAffineTransformation::AlignPrincipalAxes( *referenceImage, *floatingImage ) );
      break;
    }

  if ( xform )
    {
    if ( writeXformNativeSpaces )
      {
      cmtk::TransformChangeToSpaceAffine toNative( *xform, *referenceImage, *floatingImage );
      cmtk::XformIO::Write( &toNative.GetTransformation(), outputXformPath );
      }
    else
      {
      cmtk::XformIO::Write( xform, outputXformPath );
      }

#ifdef CMTK_USE_SQLITE
    if ( !updateDB.empty() )
      {
      try
	{
	cmtk::ImageXformDB db( updateDB );
	db.AddImagePairXform( outputXformPath, true /*always affine*/, referenceImagePath, floatingImagePath );
	}
      catch ( const cmtk::SQLite::Exception& ex )
	{
	cmtk::StdErr << "ERROR adding transformation to database: " << ex.what() << "\n";	
	}
      }
#endif
    }
  else
    {
    cmtk::StdErr << "ERROR: at least one of the two images does not have a grid-to-physical space coordinate transformation.\n";
    throw cmtk::ExitException( 1 );
    }

  return 0;
}

#include "cmtkSafeMain"
