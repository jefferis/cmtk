/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <IO/cmtkSQLite.h>

// test SQLite database creation
int
testSQLiteNew()
{
  cmtk::SQLite db( ":memory:" );
  return 0;
}

// test SQLite open of existing file
int
testSQLiteOpen()
{
  cmtk::SQLite db( CMTK_DATADIR "/empty.sqlite", true /*readOnly*/ );
  return 0;
}

// test SQLite table creation and data insertion
int
testSQLiteCreateAndInsert()
{
  cmtk::SQLite db( ":memory:" );
  db.Exec( "create table testing ( id integer primary key, data text )" );
  db.Exec( "insert into testing values ( NULL, 'test1')" );
  db.Exec( "insert into testing values ( 2, 'test2')" );
  db.Exec( "insert into testing values ( NULL, 'test3')" );
  return 0;
}

// test SQLite database query
int
testSQLiteQuery()
{
  cmtk::SQLite db( CMTK_DATADIR "/testDB.sqlite", true /*readOnly*/ );

  cmtk::SQLite::TableType table;
  db.Query( "select * from testing", table );

  return 0;
}
