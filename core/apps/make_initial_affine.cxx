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

#include <cmtkUniformVolume.h>
#include <cmtkAffineXform.h>
#include <cmtkMatrix4x4.h>
#include <cmtkMakeInitialAffineTransformation.h>

#include <cmtkVolumeIO.h>
#include <cmtkClassStream.h>
#include <cmtkClassStreamAffineXform.h>

int
main( const int argc, const char* argv[] )
{
  const char* referenceImagePath = NULL;
  const char* floatingImagePath = NULL;
  const char* outputXformPath = NULL;
  
  bool verbose = false;

  bool centerXform = false;

  int mode = 0;

  try
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Initialize affine transformation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Compute initial affine transformation by aligning centers of mass or principal axes" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] referenceImagePath floatingImagePath outputTransformationPath" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &verbose, true, "Verbose mode." );

    cl.AddSwitch( Key( 'd', "direction-vectors" ), &mode, 0, "Alignment based on image direction vectors [default]" );
    cl.AddSwitch( Key( 'c', "centers-of-mass" ), &mode, 1, "Alignment based on centers of mass (translation only)" );
    cl.AddSwitch( Key( 'p', "principal-axes" ), &mode, 2, "Alignment based on principal axes" );

    cl.AddSwitch( Key( 'C', "center-xform" ), &centerXform, true, "Set transformation center (for rotation, scale) to center of reference image." );

    cl.Parse();

    referenceImagePath = cl.GetNext();
    floatingImagePath = cl.GetNext();
    outputXformPath = cl.GetNext();
    }
  catch ( cmtk::CommandLine::Exception ex )
    {
    cmtk::StdErr << ex << "\n";
    return 1;
    }

  cmtk::UniformVolume::SmartPtr referenceImage( cmtk::VolumeIO::ReadOriented( referenceImagePath, verbose ) );
  if ( ! referenceImage )
    {
    cmtk::StdErr << "ERROR: could not read image " << referenceImagePath << "\n";
    exit( 1 );
    }

  cmtk::UniformVolume::SmartPtr floatingImage( cmtk::VolumeIO::ReadOriented( floatingImagePath, verbose ) );
  if ( ! floatingImage )
    {
    cmtk::StdErr << "ERROR: could not read image " << floatingImagePath << "\n";
    exit( 1 );
    }
  
  cmtk::AffineXform::SmartPtr xform;
  switch ( mode )
    {
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
    cmtk::ClassStream stream( outputXformPath, cmtk::ClassStream::WRITE );
    if ( stream.IsValid() )
      {
      stream.WriteString( "reference_study", referenceImagePath );
      stream.WriteString( "floating_study", floatingImagePath );
      
      stream << *xform;      
      stream.Close();
      }
    }
  else
    {
    cmtk::StdErr << "ERROR: at least one of the two images does not have a grid-to-physical space coordinate transformation.\n";
    exit( 1 );
    }

  return 0;
}
