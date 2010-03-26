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
 * The image and transformation database has three tables that store:
 * a) images and the coordinate spaces that they live in, b) coordinate
 * transformations, and from which source space to which target space they map.
 *
 * The "images" table stores all images and identifies which space they live in:
 *\code
 *  CREATE TABLE images(id INTEGER PRIMARY KEY, space INTEGER, path TEXT);
 *\endcode
 * Each image is assigned a table-unique ID. All images that live in the same
 * coordinate space (i.e., the same acquisition of the same subject) share a
 * space ID, which is the unique image ID of the first image that was added to
 * this space.
 *
 * The "xforms" table stores the file system path and properties for each coordinate
 * transformation:
 *\code
 *  CREATE TABLE xforms(id INTEGER PRIMARY KEY, path TEXT, invertible INTEGER, spacefrom INTEGER, spaceto INTEGER);
 *\endcode
 * Each transformation is assigned a table-unique ID.
 * The field "invertible" is a flag that is set if the transformation has an explicit
 * inverse (i.e., if it is affine). All transformations can be inverted nuerically, but
 * we usually prefer explicit inverses for speed and accuracy.
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
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageXformDB_h_included_
