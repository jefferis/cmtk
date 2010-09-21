/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "System/cmtkCommandLine.h"
#include "System/cmtkConsole.h"

#include "IO/cmtkClassStream.h"
#include "IO/cmtkClassStreamAffineXform.h"

#include "Base/cmtkUniformVolume.h"
#include "IO/cmtkVolumeIO.h"

bool Verbose = false;

const char* InputFilePath = NULL;
const char* OutputFilePath = NULL;
const char* OutputXformPath = NULL;

int Axis = cmtk::AXIS_Z;
int Factor = 2;
bool Padded = false;

int
main( const int argc, const char* argv[] )
{
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Split images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Split volume image into sub-images, i.e., to separate interleaved images into passes" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose operation" );

    cl.AddSwitch( Key( 'a', "axial" ), &Axis, (int)cmtk::AXIS_Z, "Interleaved axial images" );
    cl.AddSwitch( Key( 's', "sagittal" ), &Axis, (int)cmtk::AXIS_X, "Interleaved sagittal images" );
    cl.AddSwitch( Key( 'c', "coronal" ), &Axis, (int)cmtk::AXIS_Y, "Interleaved coronal images" );

    cl.AddSwitch( Key( 'x', "interleave-x" ), &Axis, (int)cmtk::AXIS_X, "Interleaved along x axis" );
    cl.AddSwitch( Key( 'y', "interleave-y" ), &Axis, (int)cmtk::AXIS_Y, "Interleaved along y axis" );
    cl.AddSwitch( Key( 'z', "interleave-z" ), &Axis, (int)cmtk::AXIS_Z, "Interleaved along z axis" );

    cl.AddOption( Key( 'f', "factor" ), &Factor, "Interleave factor. This is the number of subimages generated." );
    cl.AddSwitch( Key( 'p', "padded" ), &Padded, true, "Padded output, i.e., fill in removed slices" );

    cl.AddOption( Key( "output-xform-path" ), &OutputXformPath, "Optional path template (fprintf-style) for output affine transformation that maps input image coordinates to each output image." );

    cl.AddParameter( &InputFilePath, "InputImage", "Input image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
    cl.AddParameter( &OutputFilePath, "OutputImagePattern", "Output image path pattern. Use '%d' to substitute subimage index." )->SetProperties( cmtk::CommandLine::PROPS_IMAGE | cmtk::CommandLine::PROPS_OUTPUT );

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( InputFilePath, Verbose ) );
  if ( ! volume || ! volume->GetData() )
    {
    cmtk::StdErr << "ERROR: Could not read image " << InputFilePath << "\n";
    return 1;
    }

  for ( int i = 0; i < Factor; ++i )
    {
    cmtk::UniformVolume::SmartPtr subvolume( Padded ? volume->GetInterleavedPaddedSubVolume( Axis, Factor, i ) : volume->GetInterleavedSubVolume( Axis, Factor, i ) );

    char path[PATH_MAX];
    if ( snprintf( path, PATH_MAX, OutputFilePath, i ) > PATH_MAX )
      {
      cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
      }
    else
      {
      cmtk::VolumeIO::Write( *subvolume, path, Verbose );
      }

    if ( OutputXformPath )
      {
      if ( snprintf( path, PATH_MAX, OutputXformPath, i ) > PATH_MAX )
	{
	cmtk::StdErr << "ERROR: output path exceeds maximum path length\n";
	}
      cmtk::AffineXform xform;
      cmtk::Types::Coordinate xlate[3] = {0,0,0};
      xlate[Axis] = -i * volume->m_Delta[Axis];
      xform.SetXlate( xlate );

      cmtk::ClassStream stream( path, cmtk::ClassStream::WRITE );
      if ( stream.IsValid() )
	{
	stream << xform;
	stream.Close();
	}
      }
    }
  
  return 0;
}

