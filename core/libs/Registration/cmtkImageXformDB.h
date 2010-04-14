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
 * \code
 *   CREATE TABLE images(id INTEGER PRIMARY KEY, space INTEGER, path TEXT);
 * \endcode
 * Each image is assigned a table-unique ID. All images that live in the same
 * coordinate space (i.e., the same acquisition of the same subject) share a
 * space ID, which is the unique image ID of the first image that was added to
 * this space.
 *
 * The "xforms" table stores the file system path and properties for each coordinate
 * transformation:
 * \code
 *   CREATE TABLE xforms(id INTEGER PRIMARY KEY, path TEXT, invertible INTEGER, level INTEGER, spacefrom INTEGER, spaceto INTEGER);
 * \endcode
 * Each transformation is assigned a table-unique ID.
 * The field "invertible" is a flag that is set if the transformation has an explicit
 * inverse (i.e., if it is affine). All transformations can be inverted nuerically, but
 * we usually prefer explicit inverses for speed and accuracy.
 * The field "level" gives the refinement level of the transformation: a transformation
 * computed from two images with no initialization has level 0. A transformation computed
 * with initialization from an existing transformation has a level of that transformation 
 * plus one.
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
   */
  void AddImage( const std::string& imagePath /**!< File system path of the new image*/,
		 const std::string& spacePath = "" /**!< File system path of an existing image that lives in the same space*/ );
  
  /** Add a transformation between two images.
   *\return True if the operation was successful, false otherwise. Failure may be due to source and target image being in the same
   * space to begin with.
   */
  bool AddXform( const std::string& xformPath, /**!< File system path of the tranformation */
		 const bool invertible, /**<! Flag: does the transformation have an explicit inverse (i.e., is it affine)? */
		 const std::string& imagePathSrc, /**!< File system path of the source image */
		 const std::string& imagePathTrg /**!< File system path of the target image */ );

  /** Add a refined transformation based on an existing transformation.
   *\return True if the operation was successful, false otherwise. Failure may be due to absence of the specified original
   *  transformation in the database.
   */
  bool AddXform( const std::string& xformPath, /**!< File system path of the new tranformation */
		 const bool invertible, /**<! Flag: does the transformation have an explicit inverse (i.e., is it affine)? */
		 const std::string& xformInitPath /** Path of the transformation that was used to initialize the computation of the new transformation. */
		 const bool initInverse = false /** Flag whether the new transformation is based on the inverse of the initial transformation, i.e., from and to space need to be switched. */ );

  /// Find space that image lives in and return its key.
  Self::PrimaryKeyType FindImageSpaceID( const std::string& imagePath );

  /** Find transformation between two images.
   * Only one transformation is returned, even if more than one transformation
   * connects the two spaces.
   *
   * Forward non-invertible (i.e., nonrigid) transformations are preferred, 
   * followed by forward explicitly invertible (i.e., affine) transformations,
   * then inverses of nonrigid transformations.
   * and finally inverses of affine transformations, 
   *
   *\return True if transformation exists. If false, the two given images may still be connected via a chain of
   * multiple, concatenated transformations.
   */
  bool FindXform( const std::string& imagePathSrc, /**!< File system path of the source image */
		  const std::string& imagePathTrg, /**!< File system path of the target image */
		  std::string& xformPath, /**!< File system path of the transformation. Only valid if function returns "true." Path can be empty if both images are already in the same space. */
		  bool& inverse /**!< If this is set, the given transformation needs to be inverted. */);

  /** Get the refinement level of a transformation in the database.
   *\return The level of the given transformation: 0 for an original transformation,
   * positive for refined transformation, or -1 is transformation is not in the
   * database.
   */
  int FindXformLevel( const std::string& xformPath ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageXformDB_h_included_
