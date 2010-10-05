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

#ifndef __cmtkVolumeIO_h_included_
#define __cmtkVolumeIO_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkAnatomicalOrientation.h>

#include <IO/cmtkFileFormat.h>
#include <IO/cmtkTypedStream.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Class for input/output of 3-D image data.
 * This class is an easy-to-use wrapper around all low-level image readers and
 * writers. 
 *
 * An image can be read, for example from a Nifti/Analyze file pair, simply as follows:
 * \code
 *  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::Read( "image.hdr" ) );
 * \endcode
 *
 * Note that we typically want images to be oriented consistently, so the preferred reader call would be:
 * \code
 *  cmtk::UniformVolume::SmartPtr volume( cmtk::VolumeIO::ReadOriented( "image.hdr" );
 * \endcode
 *
 * Similarly, we can write an image simply by calling
 * \code
 *  cmtk::VolumeIO::Write( volume, "image.hdr" );
 * \endcode
 *
 * The output format is determined automatically from the file name suffix. See
 * \see VolumeIO::Write
 * for more details.
 */
class VolumeIO 
{
public:
  /// This class.
  typedef VolumeIO Self;

  /// Read volume data from filesystem.
  static UniformVolume::SmartPtr Read( const char *path, const bool verbose = false );

  /// Read grid only from filesystem.
  static UniformVolume::SmartPtr ReadGrid( const char *path, const bool verbose = false );

  /// Read grid only from filesystem and bring into standard "RAS" orientation.
  static UniformVolume::SmartPtr ReadGridOriented( const char *path, const char* orientation, const bool verbose = false );

  /// Read grid only from filesystem and bring into standard "RAS" orientation.
  static UniformVolume::SmartPtr ReadGridOriented( const char *path, const bool verbose = false )
  {
    return Self::ReadGridOriented( path, AnatomicalOrientation::ORIENTATION_STANDARD, verbose );
  }

  /** Read image from filesystem and reorient to align anatomy with coordinate axes.
   *\param orientation Three-character orientation code. The image will be brought into the orientation
   * specified by this string. Default is "RAS", i.e., the returned image will be oriented so that the
   * positive x axis is aligned with the anatomical L/R (left/right) direction, the y axis is aligned 
   * with the P/A (posterior/anterior) direction, and the y axis is aligned with the I/S (inferior/superior) 
   * direction.
   */
  static UniformVolume::SmartPtr ReadOriented( const char *path, const char* orientation, const bool verbose = false );

  /** Read image from filesystem and reorient to align anatomy with coordinate axes of standard coordinate system ("RAS").
   */
  static UniformVolume::SmartPtr ReadOriented( const char *path, const bool verbose = false )
  {
    return Self::ReadOriented( path, AnatomicalOrientation::ORIENTATION_STANDARD, verbose );
  }

  /// Write volume data to filesystem.
  static void Write( const UniformVolume& volume, const FileFormatID format, const char *path, const bool verbose = false );

  /** Write volume data to filesystem with automatic format parsing.
   * The output file format is determined automatically from the output name suffix.
   * \note Note that using ".hdr" will write a deprecated Analyze 7.5 hdr/img format pair
   *  with private extensions and questionable assumptions regarding the anatomical
   *  orientation of the image. To write a NIFTI hdr/img pair that avoids these problems,
   *  use the filename suffix ".img" (or write a single-file NIFTI using the ".nii" suffix).
   */
  static void Write( const UniformVolume& volume, const char *pathAndFormat, const bool verbose = false );

  /// Set flag for writing compressed images.
  static void SetWriteCompressedOn()
  {
    Self::WriteCompressedOn = true;
  }

  /// Clear flag for writing compressed images.
  static void SetWriteCompressedOff()
  {
    Self::WriteCompressedOn = false;
  }

  /// Get flag for writing compressed images.
  static bool GetWriteCompressed()
  {
    return Self::WriteCompressedOn;
  }

private:
  /// Global setting: write compressed images.
  static bool WriteCompressedOn;

  /** Initializer class.
   * An object of this class is automatically instantiated when a program is run.
   * Its constructor takes care of initializing VolumeIO, e.g., by evaluating the
   * IGS_WRITE_UNCOMPRESSED environment variable.
   */
  class Initializer
  {
  private:
    /// Default constructor: initialize VolumeIO settings.
    Initializer();

    /// Instance of the initializer class.
    static Initializer Instance;
  };
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVolumeIO_h_included_
