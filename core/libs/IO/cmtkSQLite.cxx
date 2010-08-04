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

#include "cmtkSQLite.h"

#include "System/cmtkConsole.h"

#include <stdlib.h>

cmtk::SQLite::SQLite
( const std::string& dbPath, const bool readOnly )
  : m_Good( false ),
    m_DebugMode( false )
{
  if ( readOnly )
    {
    this->m_Good = (sqlite3_open_v2( dbPath.c_str(), &this->m_DB, SQLITE_OPEN_READONLY, NULL /*zVFS*/ ) == SQLITE_OK);
    }
  else
    {
    this->m_Good = (sqlite3_open_v2( dbPath.c_str(), &this->m_DB, SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE, NULL /*zVFS*/ ) == SQLITE_OK);
    }
}

cmtk::SQLite::~SQLite()
{
  if ( this->m_Good )
    sqlite3_close( this->m_DB );
}

void
cmtk::SQLite::Exec( const std::string& sql )
{
  if ( ! this->Good() )
    throw Self::Exception( "Attempting operation on invalid SQLite database object" );

  if ( this->m_DebugMode )
    {
    StdErr << sql << "\n";
    }

  char* err = NULL;
  if ( sqlite3_exec( this->m_DB, sql.c_str(), NULL, NULL, &err ) != SQLITE_OK )
    {
    StdErr << "Exec " << sql << "\nSQL error: " << err << "\n";
    sqlite3_free( err );
    }
}

/// Callback for SQLite: add rows to results table.
extern "C" 
int
cmtkSQLiteQueryCallback( void* pTable, int ncols, char** rowdata, char** )
{
  cmtk::SQLite::TableType* table = static_cast<cmtk::SQLite::TableType*>( pTable );

  std::vector< std::string > tableRow( ncols );
  for ( int col = 0; col < ncols; ++col )
    {
    if ( rowdata[col] )
      tableRow[col] = std::string( rowdata[col] );
    else
      tableRow[col] = std::string( "NULL" );
    }
  table->push_back( tableRow );

  return 0;
}

void
cmtk::SQLite::Query( const std::string& sql, cmtk::SQLite::TableType& table ) const
{
  if ( ! this->Good() )
    throw Self::Exception( "Attempting operation on invalid SQLite database object" );

  if ( this->m_DebugMode )
    {
    StdErr << sql << "\n";
    }

  table.resize( 0 );

  char* err = NULL;
  if ( sqlite3_exec( this->m_DB, sql.c_str(), cmtkSQLiteQueryCallback, &table, &err ) != SQLITE_OK )
    {
    StdErr << "Query " << sql << "\nSQL error: " << err << "\n";
    sqlite3_free( err );
    }
}

bool
cmtk::SQLite::TableExists( const std::string& tableName ) const
{
  Self::TableType table;
  this->Query( "SELECT name FROM SQLite_Master WHERE name='"+tableName+"'", table );

  return table.size() && table[0].size() && (table[0][0] == tableName);
}
