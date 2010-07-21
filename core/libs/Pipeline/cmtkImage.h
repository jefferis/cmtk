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

#ifndef __cmtkImage_h_included_
#define __cmtkImage_h_included_

#include <cmtkconfig.h>

#include "Pipeline/cmtkPlane.h"

#include "Base/cmtkMacros.h"
#include "Base/cmtkTypes.h"
#include "Base/cmtkTypedArray.h"
#include "Base/cmtkScalarImage.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

/** Class for 2D image data.
 */
class Image :
  /// Inherit geometry information from Plane.
  public Plane
{
public:
  /// Create new object.
  static Image* New() { return new Image; }

  /// Scalar data type.
  igsClassParameter(ScalarDataType,DataType);

  /// Return a pointer to the object holding the image data.
  TypedArray::SmartPtr GetData();

  /** Replace existing data array by new one.
   * BEWARE! The size of the new array must match precisely the dimensions
   * as specified by the inherited Dims[] array.
   */
  void SetData( TypedArray::SmartPtr& data );

  /** Copy non-pipelined scalar image.
   */
  void SetFromScalarImage( ScalarImage *const scalarImage, const bool copyPixelData );

  /** Copy non-pipelined scalar image.
   */
  void SetFromScalarImage( const ScalarImage* scalarImage );

  /// Create copy as non-pipelined scalar image.
  ScalarImage* GetScalarImage() const;

  /** Return data at a certain grid location (pixel).
   *@param x Index of requested pixel in x-direction. Valid range is 
   * [0..Dims[0]-1].
   *@param y Index of requested pixel in y-direction. Valid range is 
   * [0..Dims[1]-1].
   *@param def Value returned when there is no valid data at the queried
   * position. This parameter defaults to zero.
   *@return The value of the pixel at position (x,y) or the value of parameter
   * def if there was no valid data for this position.
   */
  double GetDataAt( const int x, const int y, const double def = 0 );

  /** Set data at a certain grid location (pixel).
   *@param x Index of requested pixel in x-direction. Valid range is 
   * [0..Dims[0]-1].
   *@param y Index of requested pixel in y-direction. Valid range is 
   * [0..Dims[1]-1].
   *@param value Value to set pixel to.
   */
  void SetDataAt( const int x, const int y, const double value = 0 );

  /** Return data at a certain grid location (index).
   *@param index Index of requested pixel. Valid range is 
   * [0..GetNumPixels()-1].
   *@param def Value returned when there is no valid data at the queried
   * position. This parameter defaults to zero.
   *@return The value of the pixel at position (x,y) or the value of parameter
   * def if there was no valid data for this position.
   */
  double GetDataAt( const int index, const double def = 0 );

  /** Set data at a certain grid location (index).
   *@param x Index of requested pixel. Valid range is [0..GetNumPixels()-1].
   *@param value Value to set pixel to.
   */
  void SetDataAt( const int index, const double value = 0 );

  /** Return data at a certain grid location (in world coordinates).
   *@param x Location of requested pixel in x-direction.
   *@param y Location of requested pixel in y-direction. 
   *@param def Value returned when there is no valid data at the specified
   * position. This parameter defaults to zero.
   *@return The value of the pixel at position (x,y) or the value of parameter
   * "def" if there was no valid data for this position.
   */
  double GetDataAt( const double x, const double y, const double def = 0 );

protected:
  /// Default constructor.
  Image();

  /// Destructor.
  virtual ~Image() {};

  /// The actual image data.
  TypedArray::SmartPtr Data;
};

//@}

} // namespace cmtk

#endif
