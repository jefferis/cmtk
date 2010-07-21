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

#include "Base/cmtkUniformVolume.h"
#include "IO/cmtkVolumeIO.h"

#ifdef CMTK_SINGLE_COMMAND_BINARY
namespace cmtk
{
namespace apps
{
namespace reorient
{
#endif
bool Verbose = false;

const char* InFileName = NULL;
const char* OutFileName = NULL;

const char* OldOrientation = NULL;
const char* NewOrientation = cmtk::AnatomicalOrientation::ORIENTATION_STANDARD;

int
main( int argc, char* argv[] )
{
  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Reorientation" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Convert between image orientations, i.e., physically re-order pixel array and adapt stored anatomical orientation information" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] new-orientation infile outfile" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'i', "input-orientation" ), &OldOrientation, "Override input orientation. This is a three-letter code, e.g., 'RAS', 'LPI', etc." );
    cl.AddOption( Key( 'o', "output-orientation" ), &NewOrientation, "Override output orientation. Default is 'RAS', or the closest match supported by the output image file format" );
    
    if ( ! cl.Parse() ) return 1;
       
    InFileName = cl.GetNext();
    OutFileName = cl.GetNext();
    }
  catch ( const cmtk::CommandLine::Exception& ex ) 
    {
    cmtk::StdErr << ex << "\n";
    return false;
    }
  
  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::Read( InFileName, Verbose ) );
  if ( ! volume ) 
    {
    cmtk::StdErr << "ERROR: could not read image " << InFileName << "\n";
    exit( 1 );
    }

  if ( OldOrientation )
    {
    volume->m_MetaInformation[cmtk::META_IMAGE_ORIENTATION] = OldOrientation;
    volume->m_MetaInformation[cmtk::META_IMAGE_ORIENTATION_ORIGINAL] = OldOrientation;
    }
  else
    {
    OldOrientation = volume->m_MetaInformation[cmtk::META_IMAGE_ORIENTATION].c_str();
    }
  
  if ( Verbose )
    {
    cmtk::StdErr << "Reorienting from '" << OldOrientation << "' to '" << NewOrientation << "'\n";
    }

  // now reorient here in case the writer function doesn't try to write original orientation
  volume = cmtk::UniformVolume::SmartPtr( volume->GetReoriented( NewOrientation ) );
  // override original orientation to force output with desired output orientation
  volume->m_MetaInformation[cmtk::META_IMAGE_ORIENTATION_ORIGINAL] = NewOrientation;
  
  if ( Verbose )
    {
    cmtk::StdErr << "Writing to file " << OutFileName << "\n";
    }
  cmtk::VolumeIO::Write( *volume, OutFileName, Verbose );
  
  return 0;
}
#ifdef CMTK_SINGLE_COMMAND_BINARY
} // namespace reorient
} // namespace apps
} // namespace cmtk
#endif

