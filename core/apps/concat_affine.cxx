/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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
#include <System/cmtkExitException.h>

#include <Base/cmtkAffineXform.h>
#include <IO/cmtkXformIO.h>
#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkClassStreamAffineXform.h>

const char* OutputName = "concat.xform";
bool AppendToOutput = false;
bool InvertOutput = false;

int 
doMain( const int argc, const char* argv[] )
{
  cmtk::AffineXform concat;
  try 
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Concatenate affine transformations" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool computes the explicit concatenation of multiple affine coordinate transformations, each of which can be optionally inverted." );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_SYNTX, "concat_affine [options] x0 [x1 ...] \n WHERE x0 ... xN is [{-i,--inverse}] affine transformation #. "
		       "(If the first transformation in the sequence is inverted, then '--inverse' must be preceded by '--', i.e., use '-- --inverse xform.path').");

    typedef cmtk::CommandLine::Key Key;
    cl.AddOption( Key( 'o', "outfile" ), &OutputName, "Output transformation." );
    cl.AddSwitch( Key( 'a', "append" ), &AppendToOutput, true, "Append to output file [default: overwrite]." );
    cl.AddSwitch( Key( 'I', "invert-output" ), &InvertOutput, true, "Invert concatenated transformation before output [default: no]." );

    cl.Parse( argc, argv );

    bool firstXform = true;
    const char* next = cl.GetNextOptional();
    while ( next ) 
      {
      const bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
      if ( inverse ) 
	next = cl.GetNext();
      
      cmtk::Xform::SmartPtr xform( cmtk::XformIO::Read( next ) );
      if ( ! xform ) 
	{
	cmtk::StdErr << "ERROR: could not read transformation from " << next << "\n";
	throw cmtk::ExitException( 1 );
	}

      cmtk::AffineXform::SmartPtr affine( cmtk::AffineXform::SmartPtr::DynamicCastFrom( xform ) );
      if ( ! affine )
	{
	cmtk::StdErr << "ERROR: transformation " << next << " is not affine.\n";
	throw cmtk::ExitException( 1 );
	}

      if ( inverse )
	{
	affine = affine->GetInverse();
	}

      concat.Concat( *affine );
      if ( firstXform )
	{
	concat.ChangeCenter( cmtk::FixedVector<3,cmtk::Types::Coordinate>::FromPointer( affine->RetCenter() ) );
	firstXform = false;
	}
      
      next = cl.GetNextOptional();
      }
    }
  catch ( const cmtk::CommandLine::Exception& e ) 
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  cmtk::ClassStreamOutput outStream;

  if ( AppendToOutput )
    outStream.Open( OutputName, cmtk::ClassStreamOutput::MODE_APPEND );
  else
    outStream.Open( OutputName, cmtk::ClassStreamOutput::MODE_WRITE );

  if ( outStream.IsValid() )
    {
    if ( InvertOutput )
      {
      try
	{
	outStream << (*concat.GetInverse());
	}
      catch ( const cmtk::AffineXform::MatrixType::SingularMatrixException& ex )
	{
	cmtk::StdErr << "ERROR: output transformation has singular matrix and cannot be inverted\n";
	throw cmtk::ExitException( 1 );
	}
      }
    else
      {
      outStream << concat;
      }
    }
  else
    {
    cmtk::StdErr << "ERROR: could not open output file " << OutputName << "\n";
    throw cmtk::ExitException( 1 );
    }

  return 0;
}

#include "cmtkSafeMain"
