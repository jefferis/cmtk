/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#ifndef __cmtkImageInfo_h_included_
#define __cmtkImageInfo_h_included_

#include <cmtkconfig.h>

#include <stdlib.h>

#include <cmtkTypes.h>
#include <cmtkVector3D.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Image and acquisition information.
 * This class is used for handling of file format descriptions. Importing 
 * images, the import filters retrieve necessary parameters from this class
 * and store information gathered from the data files back into it. Exporting
 * images, file format parameters are defined for the output filters.
 */
class ImageInfo 
{
public:
  /// Directory in which images are located.
  char *imagepath;

  /**@name Image organisation parameters. */
  //@{
  /// Image dimensions in pixels.
  int dims[3];

  /// Header offset in bytes before pixel data begins in raw data files.
  int offset;

  /// Bytes stored per pixel.
  int bytesperpixel;

  /// Flag whether sequence of bytes should be reversed for each pixel.
  bool swapbytes;

  /// Flag indicating whether MSB is used as sign bit.
  bool signbit;

  /** Specifier describing the image data's logical data type.
   * Values are interpreted as defined by TYPE_xxx constants.
   *@see TYPE_XXX
   */
  ScalarDataType datatype;

  //@}

  /**@name Data range parameters. */
  //@{

  /// Minimum value stored in image file.
  Types::DataItem minimum;

  /// Maximum value stored in image file.
  Types::DataItem maximum;

  /// Largest value that shall still appear black.
  Types::DataItem black;

  /// Smallest value to appear white.
  Types::DataItem white;

  /// Padding value.
  union {
    unsigned char int8;
    unsigned short int16;
    unsigned int int32;
  } PaddingValue;

  /// Flag for use of padding value.
  bool Padding;
  //@}

  /**@name Coordinate calibration parameters. */
  //@{
  /** Flag indicating whether custom calibration should be used.
   * For raw and pgm image files, custom calibration is always used regardless
   * of parameter.
   */
  bool custom;

  /** Distance in mm between two consecutive slices.
   * Reading a series of images, this increment may be used if no 
   * tablepositions are stored in the image files themselves.
   */
  double slicedistance;

  /// Pixel width (x-direction) in mm.
  double calibrationx;

  /// Pixel height (y-direction) in mm.
  double calibrationy;
  
  /** Original distance in mm between two consecutive slices.
   */
  double original_slicedistance;

  /// Original pixel width (x-direction) in mm as read from image file.
  double original_calibrationx;

  /// Original pixel height (y-direction) in mm as read from image file.
  double original_calibrationy;
  
  /// Tableposition of the current image.
  double tablepos;

  /// Real-world coordinate of image origin.
  Vector3D ImagePosition;

  /// Orientations of image axes.
  Vector3D ImageOrientation[2];

  /// Original tableposition of the current image as read from file.
  double original_tablepos;
  //@}

  /// Initialize pointer members.
  ImageInfo ();

  /// Free allocated data.
  ~ImageInfo ();

  /// Return image format's standard file suffix.
  static const char* ImageSuffix ( const int format );

  /** Set bytesperpixel field.
   * According to the desiered value, also the datatype field is held 
   * up-to-date. 
   */
  void SetBytesPerPixel ( const int bpp );

  /// Set the 3D position of this image's origin.
  void SetImagePosition ( const double x, const double y, const double z );

  /// Set the orientation of this image (direction of x- and y-axis).
  void SetImageOrientation ( const double xx, const double xy, const double xz, const double yx, const double yy, const double yz );
};

//@}

} // namespace cmtk

#endif // #ifdef __cmtkImageInfo_h_included_
