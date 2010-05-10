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

#include <cmtkConsole.h>
#include <cmtkCommandLine.h>

#include <cmtkAffineXform.h>
#include <cmtkXformIO.h>
#include <cmtkClassStream.h>
#include <cmtkClassStreamAffineXform.h>

bool Verbose = false;

const char* OutputName = "concat.xform";
bool AppendToOutput = false;
bool InvertOutput = false;

int 
main( const int argc, const char* argv[] )
{
  cmtk::AffineXform concat;
  bool firstXform = true;

  try 
    {
    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Concatenate affine transformations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "[options] x0 [x1 ...] \n WHERE x0 ... xN is [{-i,--inverse}] affine transformation #" );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( 'v', "verbose" ), &Verbose, true, "Verbose mode" );

    cl.AddOption( Key( 'o', "outfile" ), &OutputName, "Output transformation." );
    cl.AddSwitch( Key( 'a', "append" ), &AppendToOutput, true, "Append to output file [default: overwrite]." );
    cl.AddSwitch( Key( 'I', "invert-output" ), &InvertOutput, true, "Invert concatenated transformation before output [default: no]." );

    cl.Parse();

    const char* next = cl.GetNextOptional();
    while ( next ) 
      {
      const bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
      if ( inverse ) 
	next = cl.GetNext();
      
      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next, Verbose ) );
      if ( ! xform ) 
	{
	cmtk::StdErr << "ERROR: could not read transformation from " << next << "\n";
	exit( 1 );
	}

      cmtk::AffineXform::SmartPtr affine( cmtk::AffineXform::SmartPtr::DynamicCastFrom( xform ) );
      if ( ! affine )
	{
	cmtk::StdErr << "ERROR: transformation " << next << " is not affine.\n";
	exit( 1 );
	}

      if ( inverse )
	{
	affine = affine->GetInverse();
	}

      concat.Concat( *affine );
      if ( firstXform )
	{
	concat.ChangeCenter( cmtk::Vector3D( affine->RetCenter() ) );
	firstXform = false;
	}
      
      next = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  cmtk::ClassStream outStream;

  if ( AppendToOutput )
    outStream.Open( OutputName, cmtk::ClassStream::APPEND );
  else
    outStream.Open( OutputName, cmtk::ClassStream::WRITE );

  if ( outStream.IsValid() )
    {
    if ( InvertOutput )
      outStream << (*concat.GetInverse());
    else
      outStream << concat;
    }
  else
    {
    cmtk::StdErr << "ERROR: could not open o0utput file " << OutputName << "\n";
    exit( 1 );
    }

  return 0;
}

