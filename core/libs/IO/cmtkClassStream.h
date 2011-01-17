/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef __cmtkClassStream_h_included_
#define __cmtkClassStream_h_included_

#include <cmtkconfig.h>

#include <IO/cmtkTypedStream.h>
#include <IO/cmtkStudy.h>

#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkParametricPlane.h>
#include <Base/cmtkLandmark.h>
#include <Base/cmtkLandmarkList.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Class for writing and reading various library classes to and from disk.
 * For all relevant objects in the Base library, input and output operators (>>
 * and <<) are defined , much like std::iostream does for basic data types.
 */
class ClassStream : 
  /// Inherit basic functionality from typed stream.
  public TypedStream 
{
public:
  /// Access modes for archives.
  typedef enum {
    /// Read-only access.
    READ = TYPEDSTREAM_READ,
    /// Write-only access.
    WRITE = TYPEDSTREAM_WRITE,
    /// Write-only access piped through zlib/gzip compression.
    WRITE_ZLIB = TYPEDSTREAM_WRITE_ZLIB,
    /// Open existing archive and append to it.
    APPEND = TYPEDSTREAM_APPEND
  } FileMode;

  /// Default constructor.
  ClassStream() : TypedStream() {}

  /** Open constructor.
   *\param filename Name of the archive to open.
   *\param mode Access mode, ie. read-only, write-only, etc.
   */
  ClassStream( const char *filename, const FileMode mode )
    : TypedStream( filename,  (TypedStreamMode) mode ) {}

  /** Open constructor for separate path and archive names.
   *\param dir Directory to open archive in.
   *\param archive Name of the archive to open.
   *\param mode Access mode, ie. read-only, write-only, etc.
   */
  ClassStream( const char *dir, const char *archive, const FileMode mode )
    : TypedStream( dir, archive, (TypedStreamMode) mode ) {}

  /** Open another archive without constructing a new object.
   */
  void Open( const char *filename, const FileMode mode ) {
    this->TypedStream::Open( filename, (TypedStreamMode) mode );
  }

  /** Open another archive in explicit directory.
   */
  void Open( const char *dir, const char *archive, const FileMode mode ) {
    this->TypedStream::Open( dir, archive, (TypedStreamMode) mode );
  }

  /** Write generic transformation object.
   * This function determines the virtual type of the transformation object
   * (spline or linear deformation) using a dynamic_cast. It then calls the
   * appropriate specialized output function.
   */
  ClassStream& operator << ( const WarpXform *warpXform );

  /** Write spline transformation object.
   * This function works on a reference rather than a pointer. It immediately
   * calls the pointer-based function defined above for the actual writing.
   */
  ClassStream& operator << ( const SplineWarpXform& splineWarpXform )
  { return (*this) << &splineWarpXform; }
  
  /// Read (spline or linear) warp transformation.
  ClassStream& operator >> ( WarpXform::SmartPtr& warpXform );

  /// Read (spline or linear) warp transformation.
  ClassStream& operator >> ( WarpXform*& warpXform );

  /// Actually read warp transformation object.
  ClassStream& Get ( WarpXform::SmartPtr& warpXform, const AffineXform* affineXform = NULL  );

  /// Actually read warp transformation object.
  ClassStream& Get ( WarpXform*& warpXform, const AffineXform* affineXform = NULL  );

private:
  /// Write actual warp transformation object.
  ClassStream& PutWarp( const WarpXform* warpXform  );

public:
  /** Write parametric plane object.
   */
  ClassStream& operator << ( const ParametricPlane *parametricPlane );

  /** Write parametric plane object.
   * This function works on a reference rather than a pointer. It immediately
   * calls the pointer-based function defined above for the actual writing.
   */
  ClassStream& operator << ( const ParametricPlane& parametricPlane )
  { return (*this) << &parametricPlane; }
  
  /// Read parametric plane.
  ClassStream& operator >> ( ParametricPlane*& parametricPlane );

  /** Write landmark list.
   */
  ClassStream& operator << ( const LandmarkList *landmarkList );

  /** Write landmark list.
   * This function works on a reference rather than a pointer. It immediately
   * calls the pointer-based function defined above for the actual writing.
   */
  ClassStream& operator << ( const LandmarkList& landmarkList )
  { return (*this) << &landmarkList; }
  
  /// Read landmark list.
  ClassStream& operator >> ( LandmarkList::SmartPtr& landmarkList );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkClassStream_h_included_
