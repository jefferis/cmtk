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

#include <cmtkVolumeFromSlices.h>

#include <math.h>

#include <cmtkVolume.h>
#include <cmtkUniformVolume.h>
#include <cmtkAffineXform.h>
#include <cmtkMathUtil.h>
#include <cmtkProgress.h>

#include <cmtkAnatomicalOrientation.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

void 
VolumeFromSlices::InitSequence
( const ImageInfo& image, const int slice_count )
{ 
  Progress::Begin( 0, slice_count, 1, "Volume image stacking" );

  Padding = false;

  if ( (image.calibrationx == 0) || (image.calibrationy==0) ) 
    {
    HandleError( "Image calibration is zero. Using 1 as default." );
    Spacing[0] = Spacing[1] = 1;
    } 
  else
    {
    Spacing[0] = image.calibrationx;
    Spacing[1] = image.calibrationy;
    }

  ImagePosition = image.ImagePosition;

  Dims[0] = image.dims[0];
  Dims[1] = image.dims[1];
  Dims[2] = slice_count;

  BytesPerPixel = image.bytesperpixel;
  SignBit = image.signbit;

  // Determine parameters for the axis swapping and reversal to occur during
  // volume assembly.
  this->DataSize = this->Dims[0] * this->Dims[1] * this->Dims[2];

  this->RawData = this->AllocDataArray( BytesPerPixel, DataSize );

  // Allocate array for axis sample points
  for ( unsigned int idx = 0; idx<3; ++idx )
    Points[idx] = Memory::AllocateArray<Types::Coordinate>( Dims[idx] );

  // Set sample points for uniform original x- and y-axis
  for ( unsigned int dim=0; dim<2; ++dim ) 
    {
    for ( int idx=0; idx < Dims[dim]; ++idx ) 
      {
      Points[dim][idx] = idx * Spacing[dim];
      }
    // Set size in axis direction.
    Size[dim] = (Dims[dim]-1) * Spacing[dim];
    }
}

void 
VolumeFromSlices::InitSequence
( const ScalarImage* image, const unsigned int numberOfSlices )
{ 
  Padding = false;

  Spacing[0] = image->GetPixelSize( AXIS_X );
  Spacing[1] = image->GetPixelSize( AXIS_Y );

  ImagePosition = image->GetImageOrigin();

  Dims[0] = image->GetDims( AXIS_X );
  Dims[1] = image->GetDims( AXIS_Y );
  Dims[2] = numberOfSlices;

  BytesPerPixel = image->GetPixelData()->GetItemSize();
  DataType = image->GetPixelData()->GetType();

  DataSize = Dims[0] * Dims[1] * Dims[2];

  VolumeDataArray = TypedArray::SmartPtr( image->GetPixelData()->NewTemplateArray( DataSize ) );
  
  // Allocate array for axis sample points
  for ( unsigned int idx = 0; idx<3; ++idx )
    Points[idx] = Memory::AllocateArray<Types::Coordinate>( Dims[idx] );
  
  // Set sample points for uniform original x- and y-axis
  for ( unsigned int dim=0; dim<2; ++dim ) 
    {
    for ( int idx=0; idx < Dims[dim]; ++idx ) 
      {
      Points[dim][idx] = idx * Spacing[dim];
      }
    // Set size in axis direction.
    Size[dim] = (Dims[dim]-1) * Spacing[dim];
    }
}

char* 
VolumeFromSlices::AllocDataArray
( const int bytesPerPixel, const int dataSize ) const 
{
  return Memory::AllocateArray<char>( bytesPerPixel * dataSize );
}

TypedArray* 
VolumeFromSlices::EncapDataArray ( const ScalarDataType dtype, void *const data, const int data_size ) 
  const
{
  return TypedArray::Create( dtype, data, data_size, IGS_FREE_ARRAY, Padding, &PaddingValue );
}

const char* 
VolumeFromSlices::FillPlane 
( unsigned int& plane, const void *data, const ImageInfo& image )
{
  const char *check = this->CheckImage( plane, image );
  if ( check ) return check;

  int bytesPerBlock = BytesPerPixel*Dims[0]*Dims[1];
  for ( int planeIdx = 0; planeIdx < image.dims[2]; ++planeIdx, ++plane ) 
    {
    const size_t offset = (planeIdx * Dims[0] * Dims[1]) * this->BytesPerPixel;
    memcpy( this->RawData+offset, ((char*)data)+offset, bytesPerBlock );
    
    if ( image.Padding ) 
      {
      Padding = true;
      memcpy( &PaddingValue, &image.PaddingValue, sizeof( PaddingValue) );
      }
    
    // set world coordinate of the plane just read
    Types::Coordinate slice_position = (ImagePosition - FirstImagePosition).EuclidNorm();
    slice_position = 1e-6 * ( MathUtil::Round( 1e+6 * slice_position) );
    Points[2][plane] = slice_position;
    }
  
  Progress::SetProgress( plane );
  return NULL;
}

const char* 
VolumeFromSlices::FillPlane 
( unsigned int& plane, const ScalarImage* image )
{
  char* rawDataPtr = static_cast<char*>( VolumeDataArray->GetDataPtr() );

  unsigned int bytesPerBlock = BytesPerPixel * Dims[0] * Dims[1];
  for ( unsigned int planeIdx = 0; planeIdx < image->GetNumberOfFrames(); ++planeIdx, ++plane ) 
    {
    const char *check = this->CheckImage( plane, image, planeIdx );
    if ( check ) return check;
    
    memcpy( rawDataPtr + bytesPerBlock * plane, image->GetPixelData()->GetDataPtr(), bytesPerBlock );
    
    // set world coordinate of the plane just read
    Types::Coordinate slicePosition = (ImagePosition - FirstImagePosition).EuclidNorm();
    slicePosition = 1e-6 * ( MathUtil::Round( 1e+6 * slicePosition) );
    Points[2][plane] = slicePosition;
    }
  
  return NULL;
}

UniformVolume* 
VolumeFromSlices::FinishVolume ( Types::Coordinate& sliceOffset, int& sliceDirection )
{
  Types::Coordinate *next_point = Points[2];

  sliceOffset = next_point[0];
  sliceDirection = MathUtil::Sign(next_point[1]-sliceOffset);
  
  Types::Coordinate previous_plane = sliceOffset;
  
  // normalize z-coordinates so that they start with zero and increase with
  // growing z-index.
  *next_point = 0;
  int idx;
  for ( idx=1, ++next_point; idx < Dims[2]; ++idx, ++next_point ) 
    {
    Types::Coordinate next_plane = *next_point;
    (*next_point) = *(next_point-1)+fabs(next_plane-previous_plane);
    previous_plane = next_plane;
    }
  
  Size[2] = *(next_point-1);
  
  // Encapsulate raw volume data.
  
  if ( !VolumeDataArray )
    VolumeDataArray = TypedArray::SmartPtr( this->EncapDataArray( SelectDataTypeInteger( BytesPerPixel, SignBit ), RawData, DataSize ) );
  
  const Types::Coordinate* aux[] = { Points[0], Points[1], Points[2] };
  UniformVolume* Result = this->ConstructVolume( Dims, Size, aux, VolumeDataArray );

  // clear reference, since now linked by volume.
  VolumeDataArray = TypedArray::SmartPtr::Null; 

  for ( idx = 0; idx<3; ++idx )
    delete[] Points[idx];

  Result->m_MetaInformation[CMTK_META_SPACE] = Result->m_MetaInformation[CMTK_META_SPACE_ORIGINAL] = "LPS";

  // actual image directions
  const Types::Coordinate spacing[3] = { (Size[0] / (Dims[0]-1)), (Size[1] / (Dims[1]-1)), (Size[2] / (Dims[2]-1)) };

  this->ImageOrientation[0] *= spacing[0] / this->ImageOrientation[0].EuclidNorm();
  this->ImageOrientation[1] *= spacing[1] / this->ImageOrientation[1].EuclidNorm();
  this->IncrementVector *= spacing[2] / this->IncrementVector.EuclidNorm();
  
  const Types::Coordinate directions[3][3] = 
    {
      { this->ImageOrientation[0][0], this->ImageOrientation[0][1], this->ImageOrientation[0][2] },
      { this->ImageOrientation[1][0], this->ImageOrientation[1][1], this->ImageOrientation[1][2] },
      { this->IncrementVector[0], this->IncrementVector[1], this->IncrementVector[2] }
    };
  
  const Matrix3x3<Types::Coordinate> m3( directions );
  Matrix4x4<Types::Coordinate> m4( m3 );
  for ( int i = 0; i < 3; ++i )
    m4[3][i] = this->FirstImagePosition.XYZ[i];

  Result->m_IndexToPhysicalMatrix = m4;
  const std::string orientationString0 = Result->GetOrientationFromDirections();
  Result->ChangeCoordinateSpace( "RAS" );

  const std::string orientationString = Result->GetOrientationFromDirections();
  Result->m_MetaInformation[CMTK_META_SPACE_UNITS_STRING] = "mm"; // seems to be implied in DICOM
  Result->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION] = orientationString;
  Result->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION_ORIGINAL] = orientationString;

  return Result;
}

const char* 
VolumeFromSlices::CheckImage
( const int plane, const ScalarImage* image, const unsigned int frame )
{
  if ( ( static_cast<unsigned int>(Dims[0]) != image->GetDims( AXIS_X ) ) || ( static_cast<unsigned int>(Dims[1]) != image->GetDims( AXIS_Y ) ) )
    return "Image size mismatch";
  
  if ( ( fabs( image->GetPixelSize( AXIS_X ) - Spacing[0] ) > CMTK_MAX_CALIB_ERROR ) ||
       ( fabs( image->GetPixelSize( AXIS_Y ) - Spacing[1] ) > CMTK_MAX_CALIB_ERROR ) )
    return "Calibration mismatch";
  
  // not too many things can go wrong for the very first slice.
  if ( plane == 0 ) 
    {
    FirstImagePosition = ImagePosition = image->GetImageOrigin( frame );
    ImageOrientation[0] = image->GetImageDirectionX();
    ImageOrientation[1] = image->GetImageDirectionY();
    return NULL;
    }

  // check whether this slice is parallel to the previous one
  for ( unsigned int dim = 0; dim<3; ++dim ) 
    {
    if ( ( fabs( ImageOrientation[0][dim] - image->GetImageDirectionX()[dim] ) > CMTK_MAX_CALIB_ERROR ) ||
	 ( fabs( ImageOrientation[1][dim] - image->GetImageDirectionY()[dim] ) > CMTK_MAX_CALIB_ERROR ) )
      return "Non-parallel image planes";
    }
  
  // Second++ slice: Compute slice-to-slice vector
  Vector3D imageToImage = image->GetImageOrigin( frame ) - ImagePosition;
  
  if ( imageToImage.MaxNorm() < CMTK_MAX_LOCALIZE_ERROR )
    {
    StdErr.printf( "Two slices at position (%f,%f,%f)\n", (float)ImagePosition[0], (float)ImagePosition[1], (float)ImagePosition[2] );
    return "Encountered two slices in identical location.";
    }
  else
    imageToImage /= imageToImage.MaxNorm();
  
  // Check whether slice-to-slice direction is orthogonal to image
  // axes.
  const Types::Coordinate scalarX = imageToImage * ImageOrientation[0];
  const Types::Coordinate scalarY = imageToImage * ImageOrientation[1];
  if ( ( fabs( scalarX ) > CMTK_MAX_ANGLE_ERROR ) || ( fabs( scalarY ) > CMTK_MAX_ANGLE_ERROR ) )
    return "Data grid must be orthogonal.";

  // if this is the second slice, save increment vector for further tests.
  if ( plane == 1 )
    IncrementVector = imageToImage;
  // otherwise, perform these tests
  else 
    {
    // Are we still going in the same direction?
    if ( (imageToImage - IncrementVector).MaxNorm() > CMTK_MAX_LOCALIZE_ERROR )
      {
      // Nope, but why? Let's give user some more hints
      if ( ( (imageToImage * IncrementVector) > 0 ) )
	// Basically same direction, so FOV has changed
	return "Field-of-view mismatch";
      else
	// Completely different direction: We're going backwards
	return "Encountered altering slice direction.";
      }
    }
  
  // Finally, save essential information about current image.
  ImagePosition = image->GetImageOrigin( frame );
  
  return NULL;
}

const char* 
VolumeFromSlices::CheckImage ( const int plane, const ImageInfo& image ) 
{
  if ( ( Dims[0] != image.dims[0] ) || ( Dims[1] != image.dims[1] ) )
    return "Image size mismatch";
  
  if ( image.bytesperpixel != BytesPerPixel )
    return "Data type mismatch";

  if ( ( fabs( image.calibrationx - Spacing[0] ) > CMTK_MAX_CALIB_ERROR ) || ( fabs( image.calibrationy - Spacing[1] ) > CMTK_MAX_CALIB_ERROR ) )
    return "Calibration mismatch";
  
  // not too many things can go wrong for the very first slice.
  if ( plane == 0 ) 
    {
    FirstImagePosition = ImagePosition = image.ImagePosition;
    ImageOrientation[0] = image.ImageOrientation[0];
    ImageOrientation[1] = image.ImageOrientation[1];
    return NULL;
    }
  
  // check whether this slice is parallel to the previous one
  for ( int axis = 0; axis<2; ++axis )
    for ( int dim = 0; dim<3; ++dim )
      {
      if ( fabs( ImageOrientation[axis][dim] - image.ImageOrientation[axis][dim] ) > CMTK_MAX_CALIB_ERROR )
	return "Non-parallel image planes";
      }
  
  // Second++ slice: Compute slice-to-slice vector
  Vector3D imageToImage = image.ImagePosition - ImagePosition;
  
  if ( imageToImage.MaxNorm() < CMTK_MAX_LOCALIZE_ERROR )
    return "Encountered two slices in identical location.";
  else
    imageToImage /= imageToImage.MaxNorm();
  
  // Check whether slice-to-slice direction is orthogonal to image
  // axes.
  const Types::Coordinate scalarX = imageToImage * ImageOrientation[0];
  const Types::Coordinate scalarY = imageToImage * ImageOrientation[1];
  if ( ( fabs( scalarX ) > CMTK_MAX_ANGLE_ERROR ) || ( fabs( scalarY ) > CMTK_MAX_ANGLE_ERROR ) )
    return "Data grid must be orthogonal.";

  // if this is the second slice, save increment vector for further tests.
  if ( plane == 1 )
    IncrementVector = imageToImage;
  // otherwise, perform these tests
  else 
    {
    // Are we still going in the same direction?
    if ( (imageToImage - IncrementVector).MaxNorm() > CMTK_MAX_LOCALIZE_ERROR )
      {
      // Nope, but why? Let's give user some more hints
      if ( ( (imageToImage * IncrementVector) > 0 ) )
	// Basically same direction, so FOV has changed
	return "Field-of-view mismatch";
      else
	// Completely different direction: We're going backwards
	return "Encountered altering slice direction.";
      }
    }
  
  // Finally, save essential information about current image.
  ImagePosition = image.ImagePosition;
  
  return NULL;
}

UniformVolume* 
VolumeFromSlices::ConstructVolume
( const int dims[3], const Types::Coordinate size[3], const Types::Coordinate *points[3], TypedArray::SmartPtr& data ) const
{
  bool isUniform = true;
  Types::Coordinate error = 0;
  for ( unsigned int dim=0; (dim<3) && isUniform; ++dim ) 
    {
    Types::Coordinate delta = points[dim][1] - points[dim][0];
    for ( int idx=2; (idx<dims[dim]) && isUniform; ++idx ) 
      {
      if ( fabs( delta - (points[dim][idx] - points[dim][idx-1]) ) > (0.0001 * delta ) )
	isUniform = false;
      error = fabs( delta - (points[dim][idx] - points[dim][idx-1]) );
      }
    }
  
  UniformVolume *result = NULL;
  if ( !isUniform )
    {
    StdErr << "WARNING: not a uniform volume (error = " << error << ")\n";
    }
  result = new UniformVolume( dims, size, data );
  
  return result;
}

} // namespace cmtk
