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
  this->ExecNoReturn( "create table images(id int primary key, path text)" );
  this->ExecNoReturn( "create table spaces(id int primary key, path text)" );
  this->ExecNoReturn( "create table xforms(id int primary key, path text)" );

  // create relationship tables
  this->ExecNoReturn( "create table imagespace(spaceid int, imageid int)" );
  this->ExecNoReturn( "create table spacexform(xformid int, spacefromid int, spacetoid int)" );
}

cmtk::ImageXformDB::PrimaryKeyType
cmtk::ImageXformDB
::AddImage( const std::string& imagePath )
{
  const std::string sql = "insert into images values (NULL,'"+imagePath+"')";
  this->ExecNoReturn( sql );
  return -1;
}

void
cmtk::ImageXformDB
::AddImageToSpace( const cmtk::ImageXformDB::PrimaryKeyType spaceKey, const cmtk::ImageXformDB::PrimaryKeyType imageKey )
{
  std::ostringstream sql;
  sql << "insert into imagespace values (" << spaceKey << "," << imageKey << ")";
  this->ExecNoReturn( sql.str() );
}

cmtk::ImageXformDB::PrimaryKeyType 
cmtk::ImageXformDB
::FindImageSpace( const cmtk::ImageXformDB::PrimaryKeyType& imageKey )
{
  return -1;
}
