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

#include <System/cmtkCommandLine.h>
#include <System/cmtkExitException.h>

#include <Registration/cmtkImageXformDB.h>

#include <string.h>

bool DebugMode = false;

int
addImages( const int argc, const char* argv[] )
{
  try
    {
    const char* dbpath = NULL;
    const char* space = NULL;
    std::vector<std::string> images;

    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Add images to the database" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This command manually adds a new image to the database, potentially assigning it to the space of an existing image. "
		       "This is useful, for example, when multiple images were acquired in spatial alignment, or when images were created using tools that do not support database updating." );

    typedef cmtk::CommandLine::Key Key;

    cl.AddParameter( &dbpath, "Database", "Database path." );
    cl.AddParameter( &space, "Space", "Image that defines coordinate space. If this image is not yet in the database, it will be added first." );
    cl.AddParameterVector( &images, "Image", "Other images that reside in the same anatomical space (same subject, same acquisition) as the specified space." );

    cl.Parse( argc, argv );

    try
      {
      cmtk::ImageXformDB db( dbpath );

      if ( DebugMode )
	db.DebugModeOn();

      db.AddImage( space );
      for ( size_t i = 0; i < images.size(); ++i )
	db.AddImage( space, images[i] );
      }
    catch ( const cmtk::ImageXformDB::Exception& e )
      {
      cmtk::StdErr << e.what() << "\n";
      return 1;
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
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

    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "List image space contents" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This command queries the database to list all images that are in the same space as the given one." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( "sort-id" ), &sortById, true, "Sort by internal database ID" );
    
    cl.AddParameter( &dbpath, "Database", "Database path." );
    cl.AddParameter( &image, "Image", "Query image path." );
    
    cl.Parse( argc, argv );
    
    try
      {
      cmtk::ImageXformDB db( dbpath, true /*readOnly*/ );

      if ( DebugMode )
	db.DebugModeOn();

      std::vector<std::string> list = db.GetSpaceImageList( db.FindImageSpaceID( image ), sortById );

      for ( size_t i = 0; i < list.size(); ++i )
	cmtk::StdOut << list[i] << "\n";
      }
    catch ( const cmtk::ImageXformDB::Exception& e )
      {
      cmtk::StdErr << e.what() << "\n";
      return 1;
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  return 0;
}

int
getXform( const int argc, const char* argv[] )
{
  try
    {
    const char* dbpath = NULL;
    const char* rimage = NULL;
    const char* fimage = NULL;
    
    bool all = false;

    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Get transformation link between two images" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This command queries the database to find transformation links that map from a given reference to a given floating image space." );

    typedef cmtk::CommandLine::Key Key;
    cl.AddSwitch( Key( "all" ), &all, true, "Print list of all transformations." );
    
    cl.AddParameter( &dbpath, "Database", "Database path." );
    cl.AddParameter( &rimage, "RefImage", "Reference image path: this is the image FROM which to map." );
    cl.AddParameter( &fimage, "FltImage", "Floating image path. this is the image TO which to map." );
    
    cl.Parse( argc, argv );
    
    try
      {
      cmtk::ImageXformDB db( dbpath, true /*readOnly*/ );

      if ( DebugMode )
	db.DebugModeOn();

      if ( all )
	{
	const std::vector<std::string> allXformsFwd = db.FindAllXforms( rimage, fimage );
	for ( size_t i = 0; i < allXformsFwd.size(); ++i )
	  {
	  cmtk::StdOut << allXformsFwd[i] << "\n";
	  }

	const std::vector<std::string> allXformsBwd = db.FindAllXforms( fimage, rimage );
	for ( size_t i = 0; i < allXformsBwd.size(); ++i )
	  {
	  cmtk::StdOut << "--inverse " << allXformsBwd[i] << "\n";
	  }
	}
      else
	{
	std::string xform;
	bool inverse;
	
	if ( db.FindXform( rimage, fimage, xform, inverse ) )
	  {
	  if ( inverse )
	    {
	    cmtk::StdOut << "--inverse ";
	    }
	  cmtk::StdOut << xform << "\n";
	  }
	else
	  {
	  cmtk::StdErr << "ERROR: no transformation can be found that maps from " << rimage << " to " << fimage << "\n";
	  return 1;
	  }
	}
      }
    catch ( const cmtk::ImageXformDB::Exception& e )
      {
      cmtk::StdErr << e.what() << "\n";
      return 1;
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    return 1;
    }

  return 0;
}

int
doMain( const int argc, const char* argv[] )
{
  int exitCode = 0;

  try
    {
    const char* command;

    cmtk::CommandLine cl;
    cl.SetProgramInfo( cmtk::CommandLine::PRG_TITLE, "Image/transformation database maintenance and query tool" );
    cl.SetProgramInfo( cmtk::CommandLine::PRG_DESCR, "This tool modifies and queries the database of images and transformations between them." );

    typedef cmtk::CommandLine::Key Key;

    cl.AddSwitch( Key( "debug" ), &DebugMode, true, "Turn on debug mode: print all SQL database commands and queries." );

    cl.AddParameter( &command, "Command", "Database command. One of \"add_images\", \"list_space\", \"get_xform\". Use '<command> --help' for detailed help." );

    cl.Parse( argc, argv );

    // get effective argc and argv for command
    const int cargc = argc-cl.GetNextIndex()+1;
    std::vector<const char*> cargv( cargc );
    cargv[0] = command;
    for ( size_t i = 1; i < cargc; ++i )
      cargv[i] = cl.GetNext();
    
    // run commands
    if ( !strcmp( command, "add_images" ) )
      exitCode = addImages( cargc, &cargv[0] );
    else if ( ! strcmp( command, "list_space" ) )
      exitCode = listSpace( cargc, &cargv[0] );
    else if ( ! strcmp( command, "get_xform" ) )
      exitCode = getXform( cargc, &cargv[0] );
    else
      {
      cmtk::StdErr << "Unknown command: " << command << "\n";
      exitCode = 1;
      }
    }
  catch ( const cmtk::CommandLine::Exception& e )
    {
    cmtk::StdErr << e << "\n";
    throw cmtk::ExitException( 1 );
    }

  return exitCode;
}

#include "cmtkSafeMain"
