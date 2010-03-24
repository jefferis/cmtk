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

#include <cmtkImageXformDB.h>

#include <cmtkConsole.h>

#include <string>
#include <sstream>
#include <cstdlib>

cmtk::ImageXformDB
::ImageXformDB( const std::string& dbPath, const bool readOnly ) 
  : cmtk::SQLite( dbPath, readOnly )
{
  this->InitNew();
}

void
cmtk::ImageXformDB
::InitNew() 
{
  // create entity tables
  this->Exec( "CREATE TABLE images(id INTEGER PRIMARY KEY, space INTEGER, path TEXT)" );
  this->Exec( "CREATE TABLE xforms(id INTEGER PRIMARY KEY, path TEXT, invertible INTEGER)" );

  // create relationship tables
  this->Exec( "CREATE TABLE spacexform(xform INTEGER, spacefrom INTEGER, spaceto INTEGER)" );
}

void
cmtk::ImageXformDB
::AddImage( const std::string& imagePath, const std::string& spacePath )
{
  if ( (spacePath == "") )
    {
    this->Exec( "INSERT INTO images (path) VALUES ('"+imagePath+"')" );
    this->Exec( "UPDATE images SET space=(SELECT id FROM images WHERE path='"+imagePath+"') WHERE path='"+imagePath+"'" );
    }
  else
    {
    PrimaryKeyType spaceKey = this->FindImageSpaceID( spacePath );
    if ( (spaceKey == 0) )
      {
      this->Exec( "INSERT INTO images (path) VALUES ('"+spacePath+"')" );
      this->Exec( "UPDATE images SET space=(SELECT space FROM images WHERE path='"+spacePath+"') WHERE path='"+spacePath+"'" );
      spaceKey = this->FindImageSpaceID( spacePath );
      }
    
    std::ostringstream sql;
    sql << "INSERT INTO images (space,path) VALUES ( " << spaceKey << ", '" << imagePath << "')";
    this->Exec( sql.str() );
    }
}

cmtk::ImageXformDB::PrimaryKeyType 
cmtk::ImageXformDB
::FindImageSpaceID( const std::string& imagePath )
{
  if ( imagePath != "" )
    {
    const std::string sql = "SELECT space FROM images WHERE path='"+imagePath+"'";
    
    SQLite::TableType table;
    this->Query( sql, table );
    
    if ( table.size() && table[0].size() )
      return atoi( table[0][0].c_str() );
    }

  return 0;
}
