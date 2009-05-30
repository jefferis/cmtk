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

#ifndef __cmtkScalarImage_h_included_
#define __cmtkScalarImage_h_included_

#include <cmtkconfig.h>

#include <cmtkMacros.h>

#include <cmtkTypes.h>
#include <cmtkTypedArray.h>
#include <cmtkVector3D.h>
#include <cmtkRectangle.h>
#include <cmtkMatrix3x3.h>
#include <cmtkInterpolator.h>

#include <assert.h>
#include <vector>

#include <cmtkSmartPtr.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Two-dimensional image with scalar pixel values.
 */
class ScalarImage
{
  /// Image dimensions.
  igsGetSetMacro2Array(unsigned int,Dims);

  /// Number of image frames for multi-frame images.
  igsGetSetMacro(unsigned int,NumberOfFrames);

  /// Pixel data.
  igsGetSetMacro(TypedArray::SmartPtr,PixelData);

  /// Pixel spacing.
  igsGetSetMacro2Array(Types::Coordinate,PixelSize);

  /// Frame-to-frame spacing.
  igsGetSetMacro(Types::Coordinate,FrameToFrameSpacing);

public:
  /// Smart pointer to ScalarImage
  typedef SmartPointer<ScalarImage> SmartPtr;

  /// Default constructor creates empty image.
  ScalarImage();

  /// Constructor with separate x and y dimensions.
  ScalarImage( const unsigned int dimsx, const unsigned int dimsy, const unsigned int numberOfFrames = 1 );

  /** Get ROI as sub-image.
   *\param roiFrom The (x,y) pixel of the original image that will make up the
   * top left (0,0) pixel of this image.
   *\param roiTo The (x,y) pixel of the original image that will make up the
   * bottom right (dimx-1,dimy-1) pixel of this image. 
   *\note Note that the ROI given by roiFrom and roiTo is an inclusive 
   * boundary, i.e., the indexed row and column will still be part of the
   * cropped image. The cropped image dimension is therefore
   * (1+roiTo[0]-roiFrom[0],1+roiTo[0]-roiFrom[0]).
   */
  ScalarImage( const ScalarImage* other, const unsigned int* roiFrom = NULL, const unsigned int* roiTo = NULL );

  /** Get ROI as sub-image.
   */
  ScalarImage( const ScalarImage* other, const IntROI2D* roi );

  /// Virtual destructor.
  virtual ~ScalarImage() {}

  /** Interpolate sub-image with affine transformation.
   */
  virtual ScalarImage* InterpolateFrom
  ( const ScalarImage* grid, const CoordinateMatrix3x3* matrix, const cmtk::Interpolators::InterpolationEnum interpolation = cmtk::Interpolators::LINEAR )
    const;

  /// Set region of interest.
  void SetROI( const IntROI2D& roi ) 
  {
    ROI = roi;
    HasROI = true;
  }

  /// Clear region of interest.
  void UnsetROI() 
  {
    HasROI = false;
  }

  /// Return cropped copy of this image.
  virtual ScalarImage* GetCropped() const 
  { 
    if ( HasROI )
      return new ScalarImage( this, &ROI );
    else
      return new ScalarImage( this );
  }
  
  /// Create pixel data array with given data type.
  void CreatePixelData( const ScalarDataType dtype ) 
  {
    PixelData = TypedArray::SmartPtr( TypedArray::Create( dtype, Dims[0] * Dims[1] * NumberOfFrames ) );
  }

  /** Origin of image in world coordinates.
   */
  Vector3D ImageOrigin;

  /// Set image origin.
  void SetImageOrigin( const Vector3D& imageOrigin ) 
  {
    ImageOrigin = imageOrigin;
  }

  /// Get image origin of given frame (default: 0).
  Vector3D GetImageOrigin( const unsigned int frame = 0 ) const;

  /** Direction of image rows relative to ImageOrigin.
   */
  igsGetSetMacro(Vector3D,ImageDirectionX);

  /** Direction of image columns relative to ImageOrigin.
   */
  igsGetSetMacro(Vector3D,ImageDirectionY);

  /** Image position from coordinate origin along axial direction.
   * This field is only meaningful if this 2D image is part of a 3D image.
   */
  igsGetSetMacro(Types::Coordinate,ImageSlicePosition);

  /** Image tilt with respect to axial position.
   * This field is only meaningful if this is an axial 2D image that is part of
   * a 3D image.
   */
  igsGetSetMacro(Types::Coordinate,ImageTiltAngle);

  /// Get number of pixels.
  unsigned int GetNumberOfPixels() const 
  { 
    return Dims[0] * Dims[1];
  }

  /// Get pixel at 2-D index.
  Types::DataItem GetPixelAt( const unsigned int i, const unsigned int j ) const 
  {
    Types::DataItem value;
    if ( PixelData->Get( value, i + Dims[0] * j ) ) return value;
    return 0;
  }

  /** Get pixel at fractional 2-D index by bilinear interpolation.
   *\return True if value is valid; false if at least one of the four neighbor
   * pixels was invalid or outside the image.
   */
  bool GetPixelAt( Types::DataItem& value, const Types::Coordinate i, const Types::Coordinate j ) const;

  /** Get pixel at fractional 2-D index by bicubic interpolation.
   *\return True if value is valid; false if at least one of the source
   * pixels was invalid or outside the image.
   */
  bool GetPixelAtCubic( Types::DataItem& value, const Types::Coordinate i, const Types::Coordinate j ) const;

  /// Set pixel at 2-D index.
  void SetPixelAt( const unsigned int i, const unsigned int j, const Types::DataItem data ) 
  {
    PixelData->Set( data, i + Dims[0] * j );
  }

  /// Clone (duplicate) this object.
  virtual ScalarImage* Clone() const;

  /** Clone (duplicate) this object.
   * This function optionally allows to reference the original pixel data 
   * rather than duplicating it. This is useful when the cloned object is only
   * generated as an intermediate object before the application of an in-place
   * image filter. In this case, not cloning the pixel data saves significant
   * heap operations. 
   *\note Note that this function is not a const member since referencing the
   * original image data requires modifying its reference counter.
   */
  virtual ScalarImage* Clone( const bool clonePixelData );

  /// Create downsampled copy of this image.
  ScalarImage* Downsample( const int factorX, int factorY = 0, ScalarImage *const target = NULL ) const;

  /** Apply binary mask to image.
   *\param maskImage Binary mask image.
   *\param invert If this flag is set, then all pixels are set to NULL where
   * the mask image is zero. Otherwise (default) the pixels are set to NULL
   * where the mask image is non-zero.
   */
  virtual void ApplyBinaryMask( const ScalarImage* maskImage, const bool invert = false );

  /**@name Image filter operators.
   * These functions perform various filter operations on the image data. For
   * each filter, two functions are implemented: GetXXXFiltered and
   * ApplyXXXFilter. The former returns a pointer to a newly allocated pixel
   * data array that contains the filtered data. The latter replaces the 
   * current object's pixel data with the filtered data.
   *\note A complete new ScalarImage object with filtered pixel data can
   * easily and efficiently be created by calling
   * object->Clone( false )->ApplyXXXFilter().
   */
  //@{
  /** Return median-filtered image data.
   *@param range Field of view of the median operator.
   */
  TypedArray* GetMedianFiltered( const byte range ) const;

  /// Replace pixel data with median-filtered data.
  ScalarImage* ApplyMedianFilter( const byte range ) 
  {
    this->SetPixelData( TypedArray::SmartPtr( this->GetMedianFiltered( range ) ) );
    return this;
  }
  
  /** Return Gauss-filtered (smoothed) image data.
   *@param stdDev Standard deviation of the smoothing kernel in world
   * coordinates. This value is internally converted to the effective kernel
   * size in pixels.
   */
  TypedArray* GetGaussFiltered( const Types::Coordinate stdDev ) const;

  /// Replace pixel data with Gauss-filtered data.
  ScalarImage* ApplyGaussFilter( const Types::Coordinate stdDev ) 
  {
    this->SetPixelData( TypedArray::SmartPtr( this->GetGaussFiltered( stdDev ) ) );
    return this;
  }

  /** Return Sobel-filtered (edge-enhanced) image data.
   * This function implements the 2-D Sobel edge operator. In particular, it
   * computed the gradient magnitude at a certian location as the square root
   * of the squared horizontal and vertical Sobel operators.
   */
  TypedArray* GetSobel2DFiltered() const;

  /// Replace pixel data with Sobel-filtered data.
  ScalarImage* ApplySobel2DFilter() 
  {
    this->SetPixelData( TypedArray::SmartPtr( this->GetSobel2DFiltered() ) );
    return this;
  }

  /** Return Laplace-filtered (edge-enhanced) image data.
   * This function implements the 2-D Laplace edge operator.
   */
  TypedArray* GetLaplace2DFiltered() const;

  /// Replace pixel data with Laplace-filtered data.
  ScalarImage* ApplyLaplace2DFilter() 
  {
    this->SetPixelData( TypedArray::SmartPtr( this->GetLaplace2DFiltered() ) );
    return this;
  }
  
  /** Return Sobel-filtered (edge-enhanced) image data.
   * This function implements the 1-D Sobel edge operator.
   */
  TypedArray* GetSobelFiltered( const bool horizontal, const bool absolute = false ) const;
  
  /// Replace pixel data with 1-D Sobel-filtered data.
  ScalarImage* ApplySobelFilter( const bool horizontal, const bool absolute = false ) 
  {
    this->SetPixelData( TypedArray::SmartPtr( this->GetSobelFiltered( horizontal, absolute ) ) );
    return this;
  }
  //@}

  /// Mirror image horizontally and/or vertically.
  void Mirror( const bool horizontal, const bool vertical );

  /// Adjust aspect ratio.
  void AdjustAspect( const bool interpolate = false );

  /// Adjust aspect ratio.
  void AdjustToIsotropic( const Types::Coordinate pixelSize, const bool interpolate = false );

  /** Project 3D coordinate onto image plane.
   *@param v Original coordinate.
   *@param i Index of projected pixel in x direction.
   *@param j Index of projected pixel in y direction.
   */
  void ProjectPixel( const Vector3D& v, unsigned int& i, unsigned int& j ) const;
  
  /// Subtract one image from another in place.
  ScalarImage& operator-=( const ScalarImage& );

  /// Print object information.
  virtual void Print() const;

private:
  /// Cropping region of interest.
  IntROI2D ROI;

  /// Flag for validity of ROI.
  bool HasROI;

  // Return filtered image data using client-provided symmetric kernels.
  TypedArray* GetFilteredData( const std::vector<Types::DataItem>& filterX, const std::vector<Types::DataItem>& filterY ) const;

  /// Adjust aspect ratio by stretching in Y-direction.
  void AdjustAspectY( const bool interpolate = false );

  /// Adjust aspect ratio by stretching in X-direction.
  void AdjustAspectX( const bool interpolate = false );

  /// Subtract one image from another.
  friend ScalarImage* operator- ( const ScalarImage&, const ScalarImage& );
};

//@{
///
#define CMTK_SCALARIMAGE_NOCLONEDATA false

#define CMTK_SCALARIMAGE_CLONEDATA true

#define CMTK_SCALARIMAGE_HORIZONTAL true

#define CMTK_SCALARIMAGE_VERTICAL false

#define CMTK_SCALARIMAGE_ABSOLUTE true

#define CMTK_SCALARIMAGE_SIGNED false

//@}

//@}

} // namespace cmtk

#endif // #ifndef __cmtkScalarImage_h_included_
