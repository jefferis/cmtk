/*
//
//  Copyright 2010 SRI International
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

#include <cmtkCommandLine.h>
#include <cmtkImageXformDB.h>

#include <cstring>

bool DebugMode = false;

int
addImage( const int argc, const char* argv[] )
{
  try
    {
    const char* dbpath = NULL;
    const char* image = NULL;
    const char* space = NULL;

    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Add an image to the database" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This command manually adds a new image to the database, potentially assigning it to the space of an existing image. "
		       "This is useful, for example, when multiple images were acquired in spatial alignment, or when images were created using tools that do not support database updating." );

    typedef cmtk::CommandLine::Key Key;

    cl.AddParameter( &dbpath, "Database", "Database path." );
    cl.AddParameter( &image, "Image", "New Image." );
    cl.AddParameter( &space, "Space", "Existing image space (i.e., a previously added image that is aligned with the new image." );

    cl.Parse();

    try
      {
      cmtk::ImageXformDB db( dbpath );

      if ( DebugMode )
	db.DebugModeOn();

      db.AddImage( image, space );
      }
    catch ( cmtk::ImageXformDB::Exception e )
      {
      cmtk::StdErr << e.what() << "\n";
      return 1;
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  return 0;
}

int
listSpace( const int argc, const char* argv[] )
{
  try
    {
    const char* dbpath = NULL;
    const char* image = NULL;

    bool sortById = false;

    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Get image space" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This command queries the database to list all images that are in the same space as the given one." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( "sort-id" ), &sortById, true, "Sort by internal database ID" );
    
    cl.AddParameter( &dbpath, "Database", "Database path." );
    cl.AddParameter( &image, "Image", "Query image path." );
    
    cl.Parse();
    
    try
      {
      cmtk::ImageXformDB db( dbpath, true /*readOnly*/ );

      if ( DebugMode )
	db.DebugModeOn();

      std::vector<std::string> list = db.GetSpaceImageList( db.FindImageSpaceID( image ), sortById );

      for ( size_t i = 0; i < list.size(); ++i )
	cmtk::StdOut << list[i] << "\n";
      }
    catch ( cmtk::ImageXformDB::Exception e )
      {
      cmtk::StdErr << e.what() << "\n";
      return 1;
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  return 0;
}

int
getXform( const int argc, const char* argv[] )
{
  return 0;
}

int
main( const int argc, const char* argv[] )
{
  int exitCode = 0;

  try
    {
    const char* command;

    cmtk::CommandLine cl( argc, argv );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image/transformation database maintenance and query tool" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool modifies and queries the database of images and transformations between them." );

    typedef cmtk::CommandLine::Key Key;

    cl.AddSwitch( Key( "debug" ), &DebugMode, true, "Turn on debug mode: print all SQL database commands and queries." );

    cl.AddParameter( &command, "Command", "Database command. One of \"add_image\", \"list_space\", \"get_xform\". Use '<command> --help' for detailed help." );

    cl.Parse();

    // get effective argc and argv for command
    const int cargc = argc-cl.GetNextIndex()+1;
    const char** cargv = argv+cl.GetNextIndex()-1;
    
    // run commands
    if ( !strcmp( command, "add_image" ) )
      exitCode = addImage( cargc, cargv );
    else if ( ! strcmp( command, "list_space" ) )
      exitCode = listSpace( cargc, cargv );
    else if ( ! strcmp( command, "get_xform" ) )
      exitCode = getXform( cargc, cargv );
    else
      {
      cmtk::StdErr << "Unknown command: " << command << "\n";
      exitCode = 1;
      }
    }
  catch ( cmtk::CommandLine::Exception e )
    {
    cmtk::StdErr << e << "\n";
    exit( 1 );
    }

  return exitCode;
}

