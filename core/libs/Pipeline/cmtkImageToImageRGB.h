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

#ifndef __cmtkImageToImageRGB_h_included_
#define __cmtkImageToImageRGB_h_included_

#include <cmtkconfig.h>

#include <cmtkMultiFilter.h>

#include <cmtkColormap.h>
#include <cmtkImage.h>
#include <cmtkImageRGB.h>
#include <cmtkMacros.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Filter to convert image to RGB image using a color lookup table.
 */
class ImageToImageRGB : 
  /// This is a filter with multiple input.
  public MultiFilter<ImageRGB> 
{
public:
  /// Create new object.
  static ImageToImageRGB* New() { return new ImageToImageRGB; }

  /// Return virtual class name.
  virtual const char *GetClassName() const { return "ImageToImageRGB"; }

  /// Enumeration of Alpha modes.
  typedef enum {
    /// Do not use alpha channel (generate 3 byte-per-pixel RGB data).
    AlphaModeNone,
    /// Use constant alpha value (opaque).
    AlphaModeConst,
    /// Use linear alpha ramp.
    AlphaModeRamp
  } igsAlphaMode;
  
  /** This flag detemines if and how an alpha (transparancy) ramp is used.
   * If this flag is said, the attached Colormap object will generate
   * transparency information in addition to the usual RGB data. As a result,
   * and RGBA image will be generated instead of a plain RGB image.
   */
  igsClassParameter(igsAlphaMode,AlphaMode);
  
  /** Lower bound of the transparency ramp.
   * All pixels with data values less than this value will be assigned complete
   * transparency (alpha=0) in the output image.
   */
  igsClassParameter(double,AlphaRampFrom);

  /** Upper bound of the transparency ramp.
   * All pixels with data values higher than this value will be assigned 
   * complete opacity (alpha=1) in the output image.
   */
  igsClassParameter(double,AlphaRampTo);

  /// Use checkerboard pattern to fill PaddingData areas.
  igsClassParameter(bool,CheckerboxPadding);

  /// Convert image to RGB image.
  virtual void Execute();

  /// Set an input image.
  void SetInput( Image *const image );

  /// Set a colormap.
  void SetColormap( Colormap *const colormap );

protected:
  /// Default constructor.
  ImageToImageRGB();

  /** Destructor.
   * Dereference all input and output objects.
   */
  virtual ~ImageToImageRGB();

private:
  /// The input image.
  Image *m_Image;

  /// The actual colormap.
  Colormap *m_Colormap;

  /// Convenience definition of the parent class type.
  typedef MultiFilter<ImageRGB> Superclass;

  /// Overwrite padded regions with checkerboard pattern.
  template<class T> void MarkPaddingData( const unsigned int dimsx, const unsigned int dimsy, T *const rgba, const TypedArray* data ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkImageToImageRGB_h_included_
