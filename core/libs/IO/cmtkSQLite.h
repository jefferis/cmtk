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

  /// Constructor: open SQLite database.
  SQLite( const std::string& dbPath, /**!< Path to the SQLite3 database file. */
	  const bool readOnly = false /**!< If this flag is set, the database is opened read-only. If false, the database is opened for read/write, and a non-existing database will be created. */);
  
  /// Destructor: close database.
  virtual ~SQLite();

  /// Execute an SQL command with no return value.
  void ExecNoReturn( const std::string& sql );

protected:
  /// Database object.
  sqlite3 *m_DB;

  /// Initialize tables in newly created database.
  virtual void InitNew() {};
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSQLite_h_included_
