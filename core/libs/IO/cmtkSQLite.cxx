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

#include <cmtkSQLite.h>

#include <cmtkConsole.h>

#include <stdlib.h>

cmtk::SQLite::SQLite
( const std::string& dbPath, const bool readOnly )
{
  if ( readOnly )
    {
    if( sqlite3_open_v2( dbPath.c_str(), &this->m_DB, SQLITE_OPEN_READONLY, NULL /*zVFS*/ ) )
      {
      cmtk::StdErr << "Can't open database read-only: " << sqlite3_errmsg( this->m_DB ) << "\n";
      sqlite3_close( this->m_DB );
      exit(1);
      }
    }
  else
    {
    if( sqlite3_open_v2( dbPath.c_str(), &this->m_DB, SQLITE_OPEN_READWRITE, NULL /*zVFS*/ ) )
      {
      if ( sqlite3_open_v2( dbPath.c_str(), &this->m_DB, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL /*zVFS*/ ) )
	{
	cmtk::StdErr << "Can't open database for writing: " << sqlite3_errmsg( this->m_DB ) << "\n";
	sqlite3_close( this->m_DB );
	exit(1);
	}
      else
	{
	// new file
	this->InitNew();
	}
      }
    }
}

void
cmtk::SQLite::ExecNoReturn( const char* sql )
{
  char* err = NULL;
  const int rc = sqlite3_exec( this->m_DB, sql, NULL, NULL, &err );
  if ( rc )
    {
    StdErr << "SQL error: " << err << "\n";
    exit(1);
    }
}


cmtk::SQLite::~SQLite()
{
  sqlite3_close( this->m_DB );
}
