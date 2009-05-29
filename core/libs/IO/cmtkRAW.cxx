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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkRAW.h>

#include <string.h>

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

#ifdef DEBUG
#  include <stdio.h>
#endif // #ifdef DEBUG

#include <cmtkCompressedStream.h>

#include <cmtkStudyImageSetRaw.h>
#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

void
RAW::Read
( const char* filename, ImageInfo& imageInfo, StudyInfo& studyInfo, const int imageIndex ) 
{
  this->SetError( 0 );

  CompressedStream stream( filename );
  if ( ! stream.IsValid() ) 
    {
    this->SetErrorMsg( "Could not open file\n%s", filename );
    return;
    }
  
  int bytesperpixel = imageInfo.bytesperpixel;
  int dim = imageInfo.dims[0] * imageInfo.dims[1];
  
  // alloc the memory for the data
  if (! this->AllocDataPtr( dim, bytesperpixel )) 
    {
    this->SetErrorMsg( "Could not allocate memory reading image." );
    return;
    }
  
  stream.Seek(imageInfo.offset,SEEK_SET); // set offset from the beginning
  
  if ( !stream.Read( DataPtr, bytesperpixel, dim ) )
    { // error occured
    this->SetErrorMsg( "Error reading %d pixels of %d byte image data from\n%s", dim, bytesperpixel, filename );
    return;
    }
  
  // swap bytes if necessary
#ifdef WORDS_BIGENDIAN
  if ( imageInfo.swapbytes && (bytesperpixel>1) )
#else
    if ( (! imageInfo.swapbytes) && (bytesperpixel>1) )
#endif
      {
	size_t f, l;
	
	for ( f=0, l=bytesperpixel-1; f<Size; 
	      f+=bytesperpixel, l+=bytesperpixel )
	  for ( int j=0; j<bytesperpixel/2; ++j ) 
	    {
	    char d = ((char*)DataPtr)[l-j];
	    ((char*)DataPtr)[l-j] = ((char*)DataPtr)[f+j];
	    ((char*)DataPtr)[f+j] = d;
	    }
      }
  
  // no tableposition in raw images.
  imageInfo.tablepos = imageInfo.slicedistance*imageIndex;
    
  imageInfo.original_tablepos = imageInfo.tablepos;
  imageInfo.original_calibrationx = imageInfo.calibrationx;
  imageInfo.original_calibrationy = imageInfo.calibrationy;

  imageInfo.SetImagePosition( 0, 0, imageInfo.tablepos );
  imageInfo.SetImageOrientation( 1, 0, 0, 0, 1, 0 );
  
  studyInfo.Clear();
}

ScalarImage* 
RAW::Read 
( const char *path, const Study* study, const int index ) const
{
  const StudyImageSetRaw *studyRaw = 
    dynamic_cast<const StudyImageSetRaw*>( study );
  if ( ! studyRaw ) 
    {
    StdErr << "Internal error: RTTI says study object is not of type StudyImageSetRaw";
    }
  
  CompressedStream stream( path );
  if ( ! stream.IsValid() ) 
    {
    StdErr.printf( "File %s could not be opened.", path );
    return NULL;
    }
  
  unsigned int bytesPerPixel = studyRaw->GetBytesPerPixel();
  unsigned int numberOfPixels = studyRaw->GetDims( AXIS_X ) * studyRaw->GetDims( AXIS_Y );
  
  if ( ! studyRaw->GetMultiFile() ) 
    {
    numberOfPixels *= studyRaw->GetDims( AXIS_Z );
    }
  
  // alloc the memory for the data
  ScalarDataType dataType = SelectDataTypeInteger( bytesPerPixel, studyRaw->GetSigned() );
  TypedArray::SmartPtr pixelDataArray = TypedArray::SmartPtr( TypedArray::Create( dataType, numberOfPixels ) );
  
  // set offset from the beginning
  stream.Seek( studyRaw->GetHeaderLength(), SEEK_SET );
  
  if ( !stream.Read( pixelDataArray->GetDataPtr(), bytesPerPixel, numberOfPixels ) ) 
    {
    // error occured
    StdErr.printf( "Error reading %d pixels of %d byte image data from\n%s", numberOfPixels, bytesPerPixel, path );
    return NULL;
    }
  
  // swap bytes if necessary
  if ( studyRaw->GetBytesPerPixel() > 1 ) 
    {
#ifdef WORDS_BIGENDIAN
    if ( !studyRaw->GetBigEndian() )
#else
      if ( studyRaw->GetBigEndian() )
#endif
	pixelDataArray->ChangeEndianness();
    }
  
  ScalarImage* image = NULL;
  
  if ( studyRaw->GetMultiFile() ) 
    {
    image = new ScalarImage( study->GetDims( AXIS_X ), study->GetDims( AXIS_Y ) );
    } 
  else
    {
    image = new ScalarImage( study->GetDims( AXIS_X ), study->GetDims( AXIS_Y ), study->GetDims( AXIS_Z ) );
    }
  
  image->SetPixelSize( study->GetCalibration( AXIS_X ), study->GetCalibration( AXIS_Y ) );
  image->SetFrameToFrameSpacing( study->GetCalibration( AXIS_Z ) );
  Types::Coordinate slicePosition = index*study->GetCalibration( AXIS_Z );
  
  image->SetImageSlicePosition( slicePosition );
  image->SetImageOrigin( Vector3D(0, 0, slicePosition ) );
  image->SetImageDirectionX( Vector3D( 1, 0, 0 ) );
  image->SetImageDirectionY( Vector3D( 0, 1, 0 ) );
  image->SetPixelData( pixelDataArray );
  
  return image;
}

void
RAW::Write
( const char* filename, const ImageInfo& imageInfo, const StudyInfo&, const int )
{
  this->SetError( 0 );

  FILE *fp;
  if ( ( fp = fopen( filename, "wb" ) ) == NULL ) 
    {
    this->SetErrorMsg( "Error opening output file." );
    return;
    }
  
  int bytesperpixel = imageInfo.bytesperpixel;
  unsigned dim = imageInfo.dims[0] * imageInfo.dims[1];
  
  int Result = 0;
  // swap bytes if necessary
#ifdef WORDS_BIGENDIAN
  if ( imageInfo.swapbytes && (bytesperpixel>1) ) {
#else
  if ( ! imageInfo.swapbytes && (bytesperpixel>1) ) 
    {
#endif
    char output[16];
    for ( size_t i=0; i<Size; i+=bytesperpixel ) 
      {
      memcpy( output, ((char*)DataPtr)+i, bytesperpixel );
      for ( int j=0; j<bytesperpixel/2; ++j ) 
	{
	char d = output[bytesperpixel-1-j];
	output[bytesperpixel-1-j] = output[j];
	output[j] = d;
	}
      Result = (1 == fwrite( output, bytesperpixel, 1, fp ));
      }
    } 
  else
    {
    Result = (dim == fwrite( DataPtr, bytesperpixel, dim, fp ));
    }
  
  fclose(fp);
  
  if ( ! Result )
    this->SetErrorMsg( "Error writing raw file." );
}

} // namespace cmtk
