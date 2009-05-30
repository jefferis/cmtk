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

#ifndef __cmtkImageRGB_h_included_
#define __cmtkImageRGB_h_included_

#include <cmtkconfig.h>

#include <cmtkPlane.h>

#include <cmtkImage.h>
#include <cmtkRGB.h>

#include <cmtkTypes.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/// Type definition for the status of the alpha-channel presence flag.
typedef enum {
  /// Image has no alpha channel (3 bytes per pixel).
  IMAGE_RGB,
  /// Image has an alpha channel (4 bytes per pixel).
  IMAGE_RGBA
} ImageAlphaToggle;

/** Class to represent ready-to-display RGB image data.
 */
class ImageRGB : 
  /// Inherit geometry from Plane.
  public Plane 
{
public:
  /// Construct new class instance.
  static ImageRGB* New();

  /// Return virtual class name.
  virtual const char *GetClassName() const { return "ImageRGB"; }

  /** Get pointer to RGB data.
   * This function checks whether a data array of the appropriate size exists.
   * If not, the old array is freed and a new one with the correct size is
   * created.
   *@param forceAlloc If this flag is true, then a data array of
   * appropriate size is allocated if it had not been done before.
   */
  const byte *GetDataPtr() const { return Data; }

  /** Get pointer to RGB data.
   * This function checks whether a data array of the appropriate size exists.
   * If not, the old array is freed and a new one with the correct size is
   * created.
   *@param forceAlloc If this flag is true, then a data array of
   * appropriate size is allocated if it had not been done before.
   */
  byte *GetDataPtr( const bool forceAlloc );

  /** Return RGBA pixel data.
   * If this is actually an RGB image only, the alpha component of the returned
   * pixel will be set to 255 (opaque).
   */
  void GetPixel( RGBA& rgb, const int index );

  /** Set RGBA pixel data.
   * If this is actually an RGB image only, the target pixel's alpha value will
   * be set to 255 (opaque).
   */
  void SetPixel( const int index, const RGBA& rgb );

  /** Return RGB pixel data.
   * If this is actually an RGBA image, the requested pixel's alpha component
   * will be ignored.
   */
  void GetPixel( RGB& rgb, const int index );

  /** Set RGB pixel data.
   * If this is actually an RGBA image, the target pixel's alpha value will be
   * set to 255 (opaque).
   */
  void SetPixel( const int index, const RGB& rgb );

  /** Set alpha channel flag.
   * This function is used to toggle the image between RGB and RGB+Alpha modes.
   */
  void SetAlphaChannel( const ImageAlphaToggle alphaChannel,
			const bool convertData = false );

  /** Return current image mode.
   */
  ImageAlphaToggle GetAlphaChannel() const { return AlphaChannel; }

  /// Return true if this image is actually grey valued.
  bool IsGreyscale() const;

protected:
  /// Default costructor.
  ImageRGB();
  
  /** Destructor.
   * Free image data array if one has been allocated.
   */
  ~ImageRGB();

private:
  /** Pointer to the RGB image data.
   * Every pixel is stored as three subsequent values R, G, B, and possibly
   * Alpha. These are all in the range 0 (black) to 255 (maximum intensity).
   */
  byte *Data;

  /** The current image mode (RGB or RGB+Alpha).
   */
  ImageAlphaToggle AlphaChannel;

  /** The number of bytes per pixel associated with the current image mode.
   *@see #AlphaChannel
   */
  unsigned int BytesPerPixel;

  /** The number of bytes allocated for the currently allocated Data array.
   * Note that this is NOT the number of pixels which, depending on the state
   * of "AlphaChannel", is only 1/4 or 1/3 of this value.
   *@see #AlphaChannel
   */
  unsigned int DataSize;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageRGB_h_included_
