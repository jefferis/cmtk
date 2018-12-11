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

#ifndef __cmtkSQLite_h_included_
#define __cmtkSQLite_h_included_

#include <cmtkconfig.h>

#include <sqlite3.h>

#include <string>
#include <vector>

#include <System/cmtkException.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Wrapper class for SQLite database.
 */
class SQLite
{
public:
  /// This class.
  typedef SQLite Self;

  /// Primary key type for the underlying database. This is used to uniquely identify table entries.
  typedef sqlite3_uint64 PrimaryKeyType;

  /// Primary key value when object is not found: this should be guaranteed to never be used by the database as an actual primary key.
  static const PrimaryKeyType NOTFOUND = static_cast<PrimaryKeyType>( -1 );

  /// Table type: matrix of strings.
  typedef std::vector< std::vector< std::string > > TableType;

  /// Exception class for class-specific error reporting.
  class Exception : 
    /// Inherit from library-level exception
    public cmtk::Exception 
  {
  public:
    /// Constructor with error message.
    Exception( const std::string& error ) : cmtk::Exception( error ) {};
  };

  /// Constructor: open SQLite database.
  SQLite( const std::string& dbPath, /*!< Path to the SQLite3 database file. */
	  const bool readOnly = false /*!< If this flag is set, the database is opened read-only. If false, the database is opened for read/write, and a non-existing database will be created. */);
  
  /// Destructor: close database.
  virtual ~SQLite();

  /// Test "good" flag.
  bool Good() const
  {
    return this->m_Good;
  }

  /// Execute an SQL command with no return value.
  void Exec( const std::string& sql );

  /// Query database and return table.
  void Query( const std::string& sql, Self::TableType& table ) const;

  /// Check if table exists.
  bool TableExists( const std::string& tableName ) const;

  /// Turn on debug mode.
  void DebugModeOn() { this->m_DebugMode = true; }

  /// Turn off debug mode.
  void DebugModeOff() { this->m_DebugMode = false; }

protected:
  /// Database object.
  mutable sqlite3 *m_DB;

  /// Flag for "good" database object.
  bool m_Good;

  /// Debug mode flag: if this is set, all executed SQL queries will be printed to standard error.
  bool m_DebugMode;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSQLite_h_included_
