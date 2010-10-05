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

/**@name Build Volume from slice images. 
 *@author Torsten Rohlfing
*/

#ifndef __cmtkVolumeFromSlices_h_included_
#define __cmtkVolumeFromSlices_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkTypes.h>
#include <Base/cmtkVolume.h>
#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkScalarImage.h>

#include <stdio.h>

/**@name Error bounds for various floating point situations.
 */
//@{
/// Maximum calibration error in mm. Provides tolerance for fp rounding.
#define CMTK_MAX_CALIB_ERROR 1e-5

/** Maximum angular error. 
 * This is the maximum difference of grid angle cosines from 90 degrees.
 * Must be more tolerant than MAX_CALIB_ERROR as the imaging devices have
 * a ROTTEN floating-point accuracy.
 */
#define CMTK_MAX_ANGLE_ERROR 1e-3

/** Maximum error allowed for image localization. 
 * Must be more tolerant than MAX_CALIB_ERROR again as the imaging devices have
 * an even WORSE than ROTTEN floating-point accuracy.
 */
#define CMTK_MAX_LOCALIZE_ERROR 1e-2
//@}

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Class for building 3D fields from slice image data.
class VolumeFromSlices 
{
public:
  /// Default constructor.
  VolumeFromSlices() : VolumeDataArray( NULL ) {}

  /// Virtual dummy destructor.
  virtual ~VolumeFromSlices() {}

protected:
  /** Start creation of new volume.
   */
  void InitSequence( const ScalarImage* image, const unsigned int numberOfSlices );

  /** Allocate memory for the 3D image data. 
   */
  virtual char* AllocDataArray( const int bytesperpixel, const int data_size ) const;

  /** Put image data into a custom data structure.
   * By default, the image data is encapsulated into a newly created
   * TypedArray object.
   */
  virtual TypedArray::SmartPtr EncapDataArray( const ScalarDataType dtype, void *const data, const int data_size ) const;
  
  /** Copy one slice of data into field.
   * This function rearranges the bytes in the given 2D image so that after
   * all slices have been copied to the 3D array, the xy-plane is always axial
   * with respect to the patient.
   */
  const char* FillPlane ( unsigned int& plane, const ScalarImage* image );

  /** Finish volume creation and free temporary storage.
   *@param sliceOffset This reference is set to the absolute slice coordinate
   * of the original image that became the first plane in the resulting volume.
   * This can be used to write images with precisely the same absolute
   * positions later.
   *@param sliceDirection This reference is set to a flag indicating whether 
   * in the original images the slice positions increased (+1) or decreased 
   * (-1) with increasing slice index.
   *@return The newly created volume object as returned by ConstructVolume().
   */
  UniformVolume::SmartPtr FinishVolume( Types::Coordinate& sliceOffset, int& sliceDirection );

  /** Finish volume creation without additional information.
   * If the additional information returned by the previous FinishVolume(...)
   * function is not reuqired, this function may be called instead.
   */
  UniformVolume::SmartPtr FinishVolume () 
  {
    Types::Coordinate dummy_c;
    int dummy_i;
    
    return FinishVolume( dummy_c, dummy_i );
  }
  
  /** Construct Volume object.
   * This function takes the geometry and data as read from the slice images.
   * Its purpose is to compose an Volume object out of these components.
   * By default, an instance of UniformVolume will be created for
   * uniformly-spaced data, while an instance of igsRectilinearVolume is
   * created for non-uniformly spaced images. Derived classes my override this
   * function to create specialized volume classes derived from the 
   * aforementioned base classes.
   *@param Dims Dimensions of the 3D data, ie. number of voxels in x-, y-, and
   * z-direction.
   *@param Size Extents of data in [mm] in x-, y-, and z-direction.
   *@param Points Positions of the grid points (voxels) with respect to the 
   * three spatial coordinates. In case the points are uniformly spaced in all
   * three dimensions, an instance of UniformVolume is created with grid
   * spacing as defined by the uniform spacings in this array. Otherwise, an
   * instance of igsRectilinearVolume is created with precisely this array as
   * its "Points" field.
   *@see igsRectilinearVolume#Points
   *@return The newly created instance of a class derived from Volume.
   *@see Volume
   */
  virtual UniformVolume::SmartPtr ConstructVolume( const DataGrid::IndexType& Dims, const UniformVolume::CoordinateVectorType& Size, const Types::Coordinate *Points[3], TypedArray::SmartPtr& Data ) const;

  /** Check image consistency.
   * This function is used to verify that all images share the same matrix 
   * size, identical pixel calibrations, and the same primitive data type. 
   * Also, slices with zero distance and changing directions of the table
   * position are detected and reported.
   *@param plane Index of this image in the sequence.
   *@param image A reference to a structure describing the next image.
   *@return A pointer to an error message, of NULL if image was okay.
   */
  const char* CheckImage ( const int plane, const ScalarImage* image, const unsigned int frame = 0 );

  /** Handle an error condition.
   * Basically, this function is intended to notify the user of errors 
   * occurring during the volume building process, such as inconsistent images.
   * By default, all errors are simply printed to the standard error output.
   * Derived classes may override this function to provide 
   * environment-specific interaction.
   *@message A textual description of the error condition.
   */
  virtual void HandleError ( const char* message ) const 
  {
    fputs ( message, stderr );
  }
  
private:
  /** Dimensions of the 3D data.
   * This array is filled with the number of voxels in x-, y-, and z-direction.
   */
  DataGrid::IndexType Dims;

  /** Size of the 3D data.
   * This array holds the extents of the 3D data in x-, y-, and z-direction.
   * All values are in real-world coordinates, ie. [mm].
   */
  UniformVolume::CoordinateVectorType Size;

  /** Axes points of the constructed volume.
   * During assembly of the 3D data, this array is filled with the positions
   * of the grid points in all three dimensions.
   */
  Types::Coordinate* Points[3];

  /// Number of voxels.
  unsigned int DataSize;

  /// Pointer to the volume data.
  char *RawData;

  /// Volume data array.
  TypedArray::SmartPtr VolumeDataArray;

  /// Number of allocated bytes per voxel.
  int BytesPerPixel;

  /// Is the data signed?
  bool SignBit;

  /// Primitive image data type.
  ScalarDataType DataType;

  /// Pixel calibration of the slice images.
  Types::Coordinate Spacing[2];

  /// X-coordinate of image origin.
  ScalarImage::SpaceVectorType FirstImagePosition;

  /// X-coordinate of image origin.
  ScalarImage::SpaceVectorType ImagePosition;

  /// X-coordinate of image origin.
  ScalarImage::SpaceVectorType ImageOrientation[2];

  /// Coordinate increment in x-direction for every block copy operation.
  int IncX;

  /// Coordinate increment in y-direction for every block copy operation.
  int IncY;

  /// Number of continuous bytes that can be copied.
  int BlockSize;

  /** Vector between the origins of subsequent images.
   * Once two images have been read, the difference of their origins in 3D 
   * space is copied to this field. The origins of subsequent slices must then
   * be in the very same direction in order to make up a rectangular 3D
   * volume.
   */
  ScalarImage::SpaceVectorType IncrementVector;

  /** Flag for pixel padding.
   * If this flag is set, PaddingValue defines a non-data value for padded
   * pixels.
   */
  bool Padding;

  /** Padding value.
   */
  union {
    unsigned char int8;
    unsigned short int16;
    unsigned int int32;
  } PaddingValue;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkVolumeFromSlices_h_included_
