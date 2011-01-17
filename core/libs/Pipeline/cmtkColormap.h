/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkColormap_h_included_
#define __cmtkColormap_h_included_

#include <cmtkconfig.h>

#include <Pipeline/cmtkPipelineObject.h>
#include <Pipeline/cmtkRGB.h>

#include <Base/cmtkTypedArray.h>

#include <IO/cmtkStudy.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class representing a configurable (HSV) colormap.
 * Ranges for H, S, and V can be defined. For all data values in a given
 * range, colors are interpolated linearly from these ranges. The number of
 * discrete colors can be chosen by the user (default is 256). For any given
 * set of parameters the lookup table is precomputed, so the actual lookup
 * of several subsequent images can be done very efficiently. 
 */
class Colormap : 
  /// Inherit basic functions from generic pipeline object.
  public PipelineObject 
{
public:
  /// Create new Colormap object.
  static Colormap* New() { return new Colormap(); }

  /// Flag for user-defined colormaps.
  igsClassParameter(bool,HaveUserMap);

  /// Two-value array defining the Hue range of the colormap.
  igsClassParameter2Array(Types::DataItem,HueRange);

  /// Two-value array defining the Saturation range of the colormap.
  igsClassParameter2Array(cmtk::Types::DataItem,SaturationRange);

  /// Two-value array defining the Value range of the colormap.
  igsClassParameter2Array(cmtk::Types::DataItem,ValueRange);

  // Gamma correction coefficient.
  igsClassParameter(cmtk::Types::DataItem,Gamma);

  /// The number of entries in the colormap, ie. the number of discrete colors.
  igsClassParameter(int,TableEntries);

  /** Two-value array defining the range of data values to map.
   * All values below the lower bound will be mapped to the first color in the
   * table while all values above the upper bound are mapped to the final color
   * in the table.
   */
  igsClassParameter2Array(cmtk::Types::DataItem,DataRange);

  /** Reverse order of table entries.
   */
  igsClassParameter(bool,Reverse);

  /** Chose one out of five predefined colormaps.
   *\param index The index of the desired standard colormap. Valid values are
   * 0 (Greylevel), 1 (Red), 2 (Green), 3 (Blue), 4 (Rainbow), 5 (Inverse
   * Rainbow). PALETTE_XXX constants are available for convenient access.
   */
  void SetStandardColormap( const int index );

  /// NULL-terminated list of standard colormap names.
  static const char *StandardColormaps[];

  /** Apply this colormap to an image to get an RGB presentation.
   *\param outPtr Pointer to a suffiently big memory segment that will hold
   * the resulting RGB data. The data will be stored as three unsigned
   * 8-bit values per pixel, representing the red, green, and blue components
   * of that pixel.
   *\param inPtr Pointer to a TypedArray object containing the data to be 
   * converted. The primitive data type can be any of the types supported by
   * TypedArray, eg. byte, short, float etc.
   *\param generateAlpha If this flag is set, a constant alpha value will be
   * generated for each pixel, resulting in 32 bits of aRGB data per pixel,
   * rather than 24 bits of RGB data. Default value for this parameter is off.
   */
  void Apply( void *const outPtr, const TypedArray* inPtr, const bool generateAlpha = false );

  /// Set colormap parameters from Study object.
  void SetFromStudy( const Study* study );

  /// Convert HSV color to RGB.
  static void HSV2RGB( RGB& rgb, Types::DataItem H, Types::DataItem S, Types::DataItem V );

protected:
  /// Default constructor.
  Colormap();

  /** Virtual destructor.
  */
  virtual ~Colormap() {}
  
  /** Execute function.
   * Called by the Update() function inherited from Object, this function
   * computes the lookup table using the parameters specified.
   */
  virtual void Execute();

private:
  /** Color lookup table.
   * This array holds the precomputed R, G, and B color components for 
   * "TableEntries" distinct data values in the range DataRange[0] throgh
   * DataRange[1].
   */
  std::vector<RGB> LookupTable;

  /// Precomputed scaling factor for data value to table index conversion.
  Types::DataItem InvDataRangeWidth;

  /** Apply table lookup for a particular primitive data type.
   * "T" is a template parameter specifying the primitive data type to lookup
   * in the color table, eg. byte, short, float etc.
   *\param outPtr Pointer to an array holding the RGB pixel data after table
   * lookup.
   *\param inPtr Pointer to the primitive data array of type T.
   *\param count Number of values in the array pointed to by inPtr. As inPtr
   * is not a TypedArray anymore, we have to make this explicit.
   *\param paddingFlag Flag for use of padding data.
   *\param paddingData Padding value. Values equal to this are ignored if "paddingFlag" is true.
   *\see Apply
   */
  template<class T>
  void ApplyPrimitive( RGB *const outPtr, const T* inPtr, const unsigned int count, const bool paddingFlag, const T paddingData ) const;

  /** Apply table lookup with constant alpha for one primitive data type.
   * "T" is a template parameter specifying the primitive data type to lookup
   * in the color table, eg. byte, short, float etc.
   *\param outPtr Pointer to an array holding the aRGB pixel data after table
   * lookup.
   *\param inPtr Pointer to the primitive data array of type T.
   *\param count Number of values in the array pointed to by inPtr. As inPtr
   * is not a TypedArray anymore, we have to make this explicit.   
   *\param paddingFlag Flag for use of padding data.
   *\param paddingData Padding value. Values equal to this are ignored if "paddingFlag" is true.
   *\see Apply
   */
  template<class T>
  void ApplyPrimitive( RGBA *const outPtr, const T* inPtr, const unsigned int count, const bool paddingFlag, const T paddingData ) const;

  /// Label color map: is system-defined by default or can be read from file.
  SegmentationLabelMap LabelColorMap;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkColormap_h_included_
