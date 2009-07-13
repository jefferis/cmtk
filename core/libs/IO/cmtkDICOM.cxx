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

#include <cmtkDICOM.h>

#ifdef CMTK_HAVE_DCMTK

#  include <cmtkDcmTags.h>
#  include <cmtkTypedArray.h>

#  include <string.h>
#  include <stdio.h>

#  include <dcdeftag.h>
#  include <didocu.h>
#  include <diutils.h>

#  ifdef CMTK_HAVE_DCMTK_JPEG
#    include <djdecode.h>
#  endif

#  ifdef HAVE_SYS_TYPES_H
#    include <sys/types.h>
#  endif

#  ifdef HAVE_TIME_H
#    include <time.h>
#  endif

#  include <memory>

#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

void
DICOM::Read
( const char* filename, ImageInfo& imageInfo, StudyInfo& studyInfo, const int imageIndex )
{
  this->SetError( 0 );
  this->FreeDataPtr();
  
  try 
    {
#ifdef CMTK_HAVE_DCMTK_JPEG
    // register global decompression codecs
    static bool decodersRegistered = false;
    if ( ! decodersRegistered ) 
      {
      DJDecoderRegistration::registerCodecs( EDC_photometricInterpretation, EUC_default, EPC_default, 1 );
      decodersRegistered = true;
      }
#endif
    
    std::auto_ptr<DcmFileFormat> fileformat( new DcmFileFormat );
    if (!fileformat.get()) 
      {
      this->SetErrorMsg( "Could not create DICOM file format object." );
      throw(1);
      }
    
    fileformat->transferInit();
    fileformat->loadFile( filename );
    fileformat->transferEnd();
    
    //    DicomImageClass::setDebugLevel( 255 );
    
    DcmDataset *dataset = fileformat->getDataset();
    if ( !dataset ) 
      {
      this->SetErrorMsg( "File format has NULL dataset." );
      throw(1);
      }
    
    std::auto_ptr<DiDocument> document( new DiDocument( dataset, dataset->getOriginalXfer() ) );
    if ( ! document.get() || ! document->good() ) 
      {
      this->SetErrorMsg( "Could not create document representation." );
      throw(2);
      }
    
    DcmElement *delem = NULL;
    Uint16 tempUint16 = 0;
    
    if ( ( delem = document->search( DCM_Rows ) ) ) 
      {
      delem->getUint16(tempUint16);
      imageInfo.dims[1]=(int)tempUint16;
      }
    
    if ( ( delem = document->search( DCM_Columns ) ) ) 
      {
      delem->getUint16(tempUint16);
      imageInfo.dims[0]=(int)tempUint16;
      }
    
    if ( ! document->getValue( DCM_NumberOfFrames, tempUint16 ) ) 
      {
      tempUint16 = 1;
      }
    imageInfo.dims[2] = tempUint16;
    
    unsigned short bitsAllocated = 0;
    if ( ( delem = document->search( DCM_BitsAllocated ) ) ) 
      {
      delem->getUint16( bitsAllocated );
      imageInfo.SetBytesPerPixel(((int)bitsAllocated+7)/8);
      } 
    else
      {
      // No "BitsAllocated" tag; use "BitsStored" instead.
      if ( ( delem = document->search( DCM_BitsStored ) ) ) 
	{
	delem->getUint16( bitsAllocated );
	imageInfo.SetBytesPerPixel(((int)bitsAllocated+7)/8);
	}
      }
    
    // get calibration from image
    bool hasPixelSpacing = (document->getValue(DCM_PixelSpacing, imageInfo.original_calibrationx, 0) > 0);
    if ( hasPixelSpacing )
      {
      if (document->getValue(DCM_PixelSpacing, imageInfo.original_calibrationy, 1) < 2) 
	{
	imageInfo.original_calibrationx = imageInfo.original_calibrationy;
	}
      } 
    else
      imageInfo.original_calibrationx = imageInfo.original_calibrationy = -1;
    
    imageInfo.signbit = false;
    if ( document->getValue( DCM_PixelRepresentation, tempUint16 ) > 0)
      imageInfo.signbit = (tempUint16 == 1);
    
    unsigned long totalImageSizePixels = imageInfo.dims[0] * imageInfo.dims[1] * imageInfo.dims[2];
    
    double rescaleIntercept, rescaleSlope;
    bool haveRescaleIntercept, haveRescaleSlope;
    if ( ! ( haveRescaleIntercept = document->getValue( DCM_RescaleIntercept, rescaleIntercept ) ) )
      rescaleIntercept = 0;
    if ( ! ( haveRescaleSlope = document->getValue( DCM_RescaleSlope, rescaleSlope ) ) )
      rescaleSlope = 1;
    
    Uint16 paddingValue;
    if ( (dataset->findAndGetUint16( DCM_PixelPaddingValue, paddingValue )).good() ) 
      {
      imageInfo.Padding = true;
      if ( imageInfo.signbit ) 
	{
	switch ( imageInfo.bytesperpixel ) 
	  {
	  case 1:
	    imageInfo.PaddingValue.int8 = static_cast<signed char>( paddingValue );
	    break;
	  case 2:
	    imageInfo.PaddingValue.int16 = static_cast<signed short>( paddingValue );
	    break;
	  case 4:
	    imageInfo.PaddingValue.int32 = static_cast<signed int>( paddingValue );
	    break;
	  default:
	    imageInfo.Padding = false;
	  }
	} 
      else 
	{
	switch ( imageInfo.bytesperpixel ) 
	  {
	  case 1:
	    imageInfo.PaddingValue.int8 = static_cast<unsigned char>( paddingValue );
	    break;
	  case 2:
	    imageInfo.PaddingValue.int16 = static_cast<unsigned short>( paddingValue );
	    break;
	  case 4:
	    imageInfo.PaddingValue.int32 = static_cast<unsigned int>( paddingValue );
	    break;
	  default:
	    imageInfo.Padding = false;
	  }
	}
      }
    
#ifdef DCM_VariablePixelData
    delem = document->search( DCM_VariablePixelData );
#else
    delem = document->search( DCM_ACR_NEMA_2C_VariablePixelData );
#endif
    if (!delem)
      delem = document->search( DCM_PixelData );
    
    if (delem) 
      {
      if ( (delem->getTag().getEVR() == EVR_OW) || (bitsAllocated > 8) ) 
	{
	Uint16 *pdata = NULL;
	delem->getUint16Array(pdata);
	this->SetDataPtr( (char*) pdata, delem->getLength() );
	if ( haveRescaleIntercept || haveRescaleSlope ) 
	  {
	  TypedArray::SmartPtr typedArray = TypedArray::SmartPtr
	    ( TypedArray::Create( TYPE_SHORT, pdata, totalImageSizePixels, IGS_NO_FREE_ARRAY, imageInfo.Padding, &imageInfo.PaddingValue ) );
#ifdef __IGNORE__
	  Types::DataItem min, max;
	  typedArray->GetRange( min, max );
	  fprintf( stderr, "before: %f %f ", min, max );
#endif
	  typedArray->Rescale( rescaleSlope, rescaleIntercept );
#ifdef __IGNORE__
	  typedArray->GetRange( min, max );
	  fprintf( stderr, "before: %f %f ", min, max );
#endif
	  imageInfo.signbit = imageInfo.signbit || (rescaleIntercept < 0);
	  }
	} 
      else 
	{
	Uint8 *pdata = NULL;
	delem->getUint8Array(pdata);
	this->SetDataPtr( (char*) pdata, delem->getLength() );
	if ( haveRescaleIntercept || haveRescaleSlope ) 
	  {
	  TypedArray::SmartPtr typedArray
	    ( TypedArray::Create( TYPE_CHAR, pdata, totalImageSizePixels, IGS_NO_FREE_ARRAY, imageInfo.Padding, &imageInfo.PaddingValue ) );
	  typedArray->Rescale( rescaleSlope, rescaleIntercept );
	  imageInfo.signbit = imageInfo.signbit || (rescaleIntercept < 0);
	  }
	}
      delem->detachValueField();
      }
    
    // now some more manual readings...
    
    // get slice spacing from multi-slice images.
    if ( ! document->getValue( DCM_SpacingBetweenSlices, imageInfo.original_slicedistance ) )
      imageInfo.original_slicedistance = 1;
    if ( ! imageInfo.custom ) 
      {
      imageInfo.slicedistance = imageInfo.original_slicedistance;
      }
    
    // get original table position from image.
    if ( ! document->getValue( DCM_SliceLocation, imageInfo.original_tablepos ) )
#ifdef DCM_Location
      if ( ! document->getValue( DCM_Location, imageInfo.original_tablepos ) )
#else
	if ( ! document->getValue( DCM_ACR_NEMA_Location, imageInfo.original_tablepos ) )
#endif
	  imageInfo.original_tablepos = imageInfo.slicedistance * imageIndex;
    
    // Use table position to set image position as long as we don't know
    // better.
    imageInfo.SetImagePosition( 0, 0, imageInfo.original_tablepos );
    
    // get original image position from file.
    const char *image_position_s = NULL;
    if ( ! document->getValue( DCM_ImagePositionPatient, image_position_s ) ) 
      {
      // ImagePositionPatient tag not present, try ImagePosition instead
#ifdef DCM_ImagePosition
      document->getValue( DCM_ImagePosition, image_position_s );
#else
      document->getValue( DCM_ACR_NEMA_ImagePosition, image_position_s );
#endif
      }
    if ( image_position_s ) 
      {
      double x, y, z;
      if ( 3 == sscanf( image_position_s,"%lf%*c%lf%*c%lf", &x, &y, &z ) ) 
	{
	imageInfo.SetImagePosition( x, y, z );
	}
      }
    
    // get original image direction from file.
    imageInfo.SetImageOrientation( 1, 0, 0, 0, 1, 0 );
    const char *image_orientation_s = NULL;
#ifdef DCM_ImageOrientation
    if ( ! document->getValue( DCM_ImageOrientation, image_orientation_s ) )
#else
      if ( ! document->getValue( DCM_ACR_NEMA_ImageOrientation, image_orientation_s ) )
#endif
	{
	// ImageOrientation tag not present, try ImageOrientationPatient
	// instead
	document->getValue( DCM_ImageOrientationPatient, image_orientation_s );
	}
    if ( image_orientation_s ) 
      {
      double xx, xy, xz, yx, yy, yz;
      if ( 6 == sscanf( image_orientation_s, "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", &xx, &xy, &xz, &yx, &yy, &yz ) ) 
	{
	imageInfo.SetImageOrientation( xx, xy, xz, yx, yy, yz );
	}
      }
    
    if ( !imageInfo.custom ) 
      {
      imageInfo.calibrationx = imageInfo.original_calibrationx;
      imageInfo.calibrationy = imageInfo.original_calibrationy;
      imageInfo.tablepos = imageInfo.original_tablepos;
      } 
    else 
      {
      imageInfo.tablepos = imageInfo.slicedistance*imageIndex;
      imageInfo.SetImagePosition( 0, 0, imageInfo.tablepos );
      imageInfo.SetImageOrientation( 1, 0, 0, 0, 1, 0 );
      }
    
    // now look up everything that goes into StudyInfo structure.
    DcmTagKey searchKey;
    const char *tmpStr;
    for ( int tag=0; NamedDcmTagTable[tag].group; ++tag ) 
      {
      searchKey.set( NamedDcmTagTable[tag].group, NamedDcmTagTable[tag].elem );
      if ( document->getValue( searchKey, tmpStr ) )
	studyInfo.SetField( NamedDcmTagTable[tag].id, strdup( tmpStr ) );
      }
    
    if ( ! DataPtr ) 
      {
      this->SetErrorMsg( "Could not read pixel data from image file\n%s", filename );
      }
  }
  
  catch ( int where ) 
    {
    switch (where) 
      {
      case 0: 
	this->SetErrorMsg( "Could not open image file\n%s", filename );
	return;
      }
    StdErr << "ERROR: " << this->ErrorMsg << "\n";
    }
}

void
DICOM::Write
( const char* filename, const ImageInfo& imageInfo, const StudyInfo& studyInfo, const int anonymize )
{
  DcmDataset *dataset = new DcmDataset;
  DcmElement *elem;

  char tmp_str[128];
  snprintf( tmp_str, sizeof( tmp_str ), "%.6f\\%.6f", (float) imageInfo.calibrationx, (float) imageInfo.calibrationy );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0030 ) )) -> putString( tmp_str );
  dataset->insert( elem );
  
  snprintf( tmp_str, sizeof( tmp_str ), "%.6f", imageInfo.tablepos );
  (elem = newDicomElement( DcmTag( 0x0020, 0x0050 ) )) -> putString( tmp_str );
  dataset->insert( elem );
  
  // Put modifying device / date / time information
  (elem = newDicomElement( DcmTag( 0x0020, 0x3401 ) )) -> putString( StudyInfo::ModDevID() );
  dataset->insert( elem );
  
#if defined(HAVE_TIMET) && defined(HAVE_LOCALTIME) && defined(HAVE_STRFTIME)
  time_t systime;
  time( &systime );
  struct tm* systm = localtime( &systime );
  
  strftime( tmp_str, 128, "%Y.%m.%d", systm );
  (elem = newDicomElement( DcmTag( 0x0020, 0x3403 ) )) -> putString( tmp_str );
  dataset->insert( elem );

  strftime( tmp_str, 128, "%H:%M:%S.0000", systm );
  (elem = newDicomElement( DcmTag( 0x0020, 0x3405 ) )) -> putString( tmp_str );
  dataset->insert( elem );
#endif

  // Put image information

  (elem = newDicomElement( DcmTag( 0x0028, 0x0002 ) )) -> putUint16( 1 ); // Samples Per Pixel
  dataset->insert( elem );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0005 ) )) -> putUint16( 2 ); // Image Dimensions
  dataset->insert( elem );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0010 ) )) -> putUint16( imageInfo.dims[1] );
  dataset->insert( elem );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0011 ) )) -> putUint16( imageInfo.dims[0] );
  dataset->insert( elem );
  
  (elem = newDicomElement( DcmTag( 0x0028, 0x0004 ) )) -> putString( "MONOCHROME" );
  dataset->insert( elem );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0100 ) )) -> putUint16( 8*imageInfo.bytesperpixel );
  dataset->insert( elem );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0101 ) )) -> putUint16( 8*imageInfo.bytesperpixel );
  dataset->insert( elem );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0102 ) )) -> putUint16( 8*imageInfo.bytesperpixel-1 );
  dataset->insert( elem );
  (elem = newDicomElement( DcmTag( 0x0028, 0x0103 ) )) -> putUint16( 0 );
  dataset->insert( elem );

  DcmElement *delem = newDicomElement( DcmTag( 0x7fe0, 0x0010 ) );
  if (delem)
    {
    if ( imageInfo.bytesperpixel == 2 )
      delem->putUint16Array( (Uint16*)DataPtr, imageInfo.dims[0] * imageInfo.dims[1] );
    else
      delem->putUint8Array( (Uint8*) DataPtr, imageInfo.dims[0] * imageInfo.dims[1] );
    dataset->insert(delem);
    }
  
  for ( int tag=0; NamedDcmTagTable[tag].group; ++tag ) 
    {
    elem = newDicomElement( DcmTag( NamedDcmTagTable[tag].group, NamedDcmTagTable[tag].elem ) );
    const char *field_value = studyInfo.GetField( NamedDcmTagTable[tag].id, anonymize );
    
    char buffer[32];
    switch ( NamedDcmTagTable[tag].id ) 
      {
      case INFO_RELIMGPOSITION:
      case INFO_RELIMGPOSITPAT:
	snprintf( buffer, sizeof( buffer ), "0\\0\\%.6f", imageInfo.tablepos );
	elem -> putString( buffer );
	break;
      case INFO_RELIMGORIENT:
      case INFO_RELIMGORIENTPAT:
	elem -> putString( "1\\0\\0\\0\\1\\0" );
	break;
      default:
	elem -> putString( field_value );
      }
    dataset->insert( elem );
    }
  
  DcmFileFormat *fileformat = new DcmFileFormat( dataset );
  
  if (!fileformat) 
    {
    this->SetErrorMsg( "Could not create DICOM file format object." );
    }
  
  fileformat->transferInit();
  fileformat->saveFile( filename, EXS_LittleEndianImplicit, EET_ExplicitLength );
  fileformat->transferEnd();
  
  delete fileformat;
}
 
ScalarImage* 
DICOM::Read 
( const char *path, const Study* study, const int index ) const
{
  ScalarImage* image = NULL;

#ifdef CMTK_HAVE_DCMTK_JPEG
  // register global decompression codecs
  static bool decodersRegistered = false;
  if ( ! decodersRegistered ) 
    {
    DJDecoderRegistration::registerCodecs( EDC_photometricInterpretation, EUC_default, EPC_default, 1 );
    decodersRegistered = true;
    }
#endif
  
  std::auto_ptr<DcmFileFormat> fileformat( new DcmFileFormat );
  if (!fileformat.get()) 
    {
    StdErr << "ERROR: Could not create DICOM file format object.\n";
    return NULL;
    }
  
  fileformat->transferInit();
  fileformat->loadFile( path );
  fileformat->transferEnd();

  //DicomImageClass::setDebugLevel( 255 );

  DcmDataset *dataset = fileformat->getDataset();

  if ( !dataset ) 
    {
    StdErr << "ERROR: File format has NULL dataset.\n";
    return NULL;
    }
  
  std::auto_ptr<DiDocument> document( new DiDocument( dataset, dataset->getOriginalXfer(), CIF_AcrNemaCompatibility ) );
  if ( ! document.get() || ! document->good() ) 
    {
    StdErr << "ERROR: Could not create document representation.\n";
    return NULL;
    }
  
  Uint16 dimsX, dimsY;
  if ( ! ( document->getValue( DCM_Columns, dimsX ) && document->getValue( DCM_Rows, dimsY ) ) ) 
    {
    StdErr << "ERROR: File " << path << " has no DCM_Columns/DCM_Rows tags.\n";
    return NULL; // no image dimensions, nothing we can do.
    }
  
  Uint16 numberOfFrames = 1;
  if ( ! document->getValue( DCM_NumberOfFrames, numberOfFrames ) )
    numberOfFrames = 1;
  image = new ScalarImage( dimsX, dimsY, numberOfFrames );
  
  // get calibration from image
  double calibrationX = -1, calibrationY = -1;
  bool hasPixelSpacing = (document->getValue(DCM_PixelSpacing, calibrationX, 0) > 0);
  if ( hasPixelSpacing ) 
    {
    if ( document->getValue(DCM_PixelSpacing, calibrationY, 1) < 2)
      calibrationY = calibrationX;
    } 
  else
    calibrationX = calibrationY = -1;
  image->SetPixelSize( calibrationX, calibrationY );
  
  Uint16 bitsAllocated = 0;
  if ( ! document->getValue( DCM_BitsAllocated, bitsAllocated ) )
    // No "BitsAllocated" tag; use "BitsStored" instead.
    document->getValue( DCM_BitsStored, bitsAllocated );
  
  bool pixelDataSigned = false;
  Uint16 pixelRepresentation = 0;
  if ( document->getValue( DCM_PixelRepresentation, pixelRepresentation ) > 0)
    pixelDataSigned = (pixelRepresentation == 1);
  
  unsigned long totalImageSizePixels = image->GetDims( AXIS_X ) * image->GetDims( AXIS_Y ) * image->GetNumberOfFrames();

  double rescaleIntercept, rescaleSlope;
  bool haveRescaleIntercept, haveRescaleSlope;
  if ( ! ( haveRescaleIntercept = document->getValue( DCM_RescaleIntercept, rescaleIntercept ) ) )
    rescaleIntercept = 0;
  if ( ! ( haveRescaleSlope = document->getValue( DCM_RescaleSlope, rescaleSlope ) ) )
    rescaleSlope = 1;
  
  pixelDataSigned = pixelDataSigned || (rescaleIntercept < 0);
  
  TypedArray* pixelDataArray( NULL );

#ifdef DCM_VariablePixelData
  DcmElement *delem = document->search( DCM_VariablePixelData );
#else
  DcmElement *delem = document->search( DCM_ACR_NEMA_2C_VariablePixelData );
#endif
  if (!delem)
    delem = document->search( DCM_PixelData );
    
  if (delem) 
    {
    if ( (delem->getTag().getEVR() == EVR_OW) || (bitsAllocated > 8) ) 
      {
      Uint16 *pdata = NULL;
      delem->getUint16Array(pdata);
      if ( pixelDataSigned ) 
	{
	pixelDataArray = TypedArray::Create( TYPE_SHORT, pdata, totalImageSizePixels );
	} 
      else
	{
	pixelDataArray = TypedArray::Create( TYPE_USHORT, pdata, totalImageSizePixels );
	}
      } 
    else
      {
      Uint8 *pdata = NULL;
      delem->getUint8Array(pdata);
      if ( pixelDataSigned ) 
	{
	pixelDataArray = TypedArray::Create( TYPE_CHAR, pdata, totalImageSizePixels );
	} 
      else 
	{
	pixelDataArray = TypedArray::Create( TYPE_BYTE, pdata, totalImageSizePixels );
	}
      }
    delem->detachValueField();
    }
  
  if ( ! pixelDataArray ) 
    {
    StdErr.printf( "Could not read pixel data from image file\n%s", path );
    }
  
  Uint16 paddingValue = 0;
  if ( ( dataset->findAndGetUint16( DCM_PixelPaddingValue, paddingValue )).good() )
    pixelDataArray->SetPaddingValue( paddingValue );
  
  if ( haveRescaleIntercept || haveRescaleSlope )
    pixelDataArray->Rescale( rescaleSlope, rescaleIntercept );

  image->SetPixelData( TypedArray::SmartPtr( pixelDataArray ) );

  // now some more manual readings...

    // get slice spacing from multi-slice images.
  double frameToFrame = 0;
  if ( document->getValue( DCM_SpacingBetweenSlices, frameToFrame ) )
    image->SetFrameToFrameSpacing( frameToFrame );

    // get original table position from image.
  double sliceLocation = 0;
  if ( ! document->getValue( DCM_SliceLocation, sliceLocation ) ) 
    {
#ifdef DCM_Location
    document->getValue( DCM_Location, sliceLocation );
#else
    document->getValue( DCM_ACR_NEMA_Location, sliceLocation );
#endif
    }
  image->SetImageSlicePosition( sliceLocation );
  
  // Use table position to set image position as long as we don't know
  // better.
  Vector3D imageOrigin( 0.0, 0.0, sliceLocation );
      
  // get original image position from file.
  const char *image_position_s = NULL;
  if ( ! document->getValue( DCM_ImagePositionPatient, image_position_s ) ) 
    {
    // ImagePositionPatient tag not present, try ImagePosition instead
#ifdef DCM_ImagePosition
    document->getValue( DCM_ImagePosition, image_position_s );
#else
    document->getValue( DCM_ACR_NEMA_ImagePosition, image_position_s );
#endif
    }
  if ( image_position_s ) 
    {
    double x, y, z;
    if ( 3 == sscanf( image_position_s,"%lf%*c%lf%*c%lf", &x, &y, &z ) ) 
      {
      imageOrigin.Set( x, y, z );
      }
    }
  
  image->SetImageOrigin( imageOrigin );
  
  // get original image direction from file.
  Vector3D imageDirectionX( 1, 0, 0 );
  Vector3D imageDirectionY( 0, 1, 0 );
  const char *image_orientation_s = NULL;
  document->getValue( DCM_ImageOrientationPatient, image_orientation_s );
  if ( image_orientation_s ) 
    {
    double xx, xy, xz, yx, yy, yz;
    if ( 6 == sscanf( image_orientation_s, "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", &xx, &xy, &xz, &yx, &yy, &yz ) ) 
      {
      imageDirectionX.Set( xx, xy, xz );
      imageDirectionY.Set( yx, yy, yz );
      }
    }
  
  image->SetImageDirectionX( imageDirectionX );
  image->SetImageDirectionY( imageDirectionY );

  if ( study && study->GetCustomCalibration() ) 
    {
    image->SetPixelSize( study->GetCalibration( AXIS_X ), study->GetCalibration( AXIS_Y ) );
    
    Types::Coordinate slicePosition = index*study->GetCalibration( AXIS_Z );
    image->SetImageSlicePosition( slicePosition );
    image->SetImageOrigin( Vector3D(0, 0, slicePosition ) );
    
    image->SetImageDirectionX( Vector3D( 1, 0, 0 ) );
    image->SetImageDirectionY( Vector3D( 0, 1, 0 ) );
    }

  return image;
}

//@}

} // namespace cmtk

#else // #ifdef CMTK_HAVE_DCMTK
namespace
cmtk
{

/** \addtogroup IO */
//@{

void DICOM::Read ( const char*, ImageInfo&, StudyInfo&, const int ) {}

ScalarImage* 
DICOM::Read( const char*, const Study*, const int ) const
{ return NULL; }

void DICOM::Write ( const char*, const ImageInfo&, const StudyInfo&, const int ) {}

//@}

} // namespace cmtk

#endif // #ifdef CMTK_HAVE_DCMTK
