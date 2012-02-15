/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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
#include <Base/cmtkMatrix4x4.h>

#include <IO/cmtkXformIO.h>

int
doMain( const int argc, const char* argv[] )
{
  const char* inputFileName = NULL;
  
  bool transpose = false;
  
  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Degrees of freedom to matrix" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "Convert affine transformation from degrees-of-freedom representation to matrix form" );

    cl.AddParameter( &inputFileName, "InputPath", "Path to input affine transformation." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );  

    typedef cmtk::CommandLine::Key Key;
    cl.BeginGroup( "Output", "Output Options" );
    cl.AddSwitch( Key( "transpose" ), &transpose, true, "Print transpose of transformation matrix." );    
    cl.EndGroup();

    cl.Parse( argc, argv );
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e;
    throw cmtk::ExitException( 1 );
    }

  if ( inputFileName )
    {
    cmtk::AffineXform::SmartConstPtr affineXform = cmtk::AffineXform::SmartPtr::DynamicCastFrom( cmtk::XformIO::Read( inputFileName ) );

    if ( affineXform )
      {
      for ( size_t j = 0; j < 4; ++j ) 
	{
	for ( size_t i = 0; i < 4; ++i ) 
	  {
	  if ( transpose )
	    {
	    cmtk::StdOut << affineXform->Matrix[i][j] << "\t";
	    }
	  else
	    {
	    cmtk::StdOut << affineXform->Matrix[j][i] << "\t";
	    }
	  }
	cmtk::StdOut << "\n";
	}
      }
    }

  return 0;
}

#include "cmtkSafeMain"
