/*
//
//  Copyright 2012 SRI International
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
//  $Revision: 1652 $
//
//  $LastChangedDate: 2010-05-14 14:45:52 -0700 (Fri, 14 May 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <cmtkconfig.h>

#include <System/cmtkCommandLine.h>
#include <System/cmtkConsole.h>

#include <string>

int
doMain( const int argc, const char* argv[] )
{

  std::string programTitle;
  std::string programDescription;

  try
    {
    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, programTitle );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, programDescription );

    typedef cmtk::CommandLine::Key Key;
//    cl.AddOption( Key( "inversion-tolerance" ), &inversionTolerance, "Numerical tolerance of B-spline inversion in mm. Smaller values will lead to more accurate inversion, but may increase failure rate." );
//    cl.AddOption( Key( "source-image" ), &sourceImagePath, "Set source image of the transformation (i.e., the image that the transformation maps points FROM) to correct for differences in orientation and coordinate space." );
//    cl.AddOption( Key( "target-image" ), &targetImagePath, "Set target image of the transformation (i.e., the image that the transformation maps points TO) to correct for differences in orientation and coordinate space." );

//    cl.AddParameterVector( &inputXformPaths, "XformList", "List of concatenated transformations. Insert '--inverse' to use the inverse of the transformation listed next." )->SetProperties( cmtk::CommandLine::PROPS_XFORM );  

    cl.Parse( argc, argv );
    }
  catch ( cmtk::CommandLine::Exception ex )
    {
    cmtk::StdErr << ex << "\n";
    throw cmtk::ExitException( 1 );
    }

  return 0;
}

#include "cmtkSafeMain"
