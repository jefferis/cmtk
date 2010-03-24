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

#ifndef __cmtkImageXformDB_h_included_
#define __cmtkImageXformDB_h_included_

#include <cmtkconfig.h>

#include <cmtkSQLite.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/** Class for image and transformation database.
 */
class ImageXformDB
/// Inherit from SQLite wrapper class.
  : public SQLite
{
public:
  /// This class.
  typedef ImageXformDB Self;

  /// Parent class.
  typedef SQLite Superclass;

  /// Constructor: open ImageXformDB database.
  ImageXformDB( const std::string& dbPath, /**!< Path to the database file. */
		const bool readOnly = false /**!< If this flag is set, the database is opened read-only. If false, the database is opened for read/write, and a non-existing database will be created. */);

  /** Add an image to a coordinate space, each identified by its file system path.
   *\return Key of newly entered image.
   */
  void AddImage( const std::string& imagePath /**!< File system path of the new image*/,
		 const std::string& spacePath = "" /**!< File system path of an existing image that lives in the same space*/ );
  
  /// Find space that image lives in and return its key.
  Self::PrimaryKeyType FindImageSpaceID( const std::string& imagePath );
  
protected:
  /// Initialize tables in newly created database.
  virtual void InitNew();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageXformDB_h_included_
