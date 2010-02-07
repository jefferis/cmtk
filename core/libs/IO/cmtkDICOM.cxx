/*
//
//  Copyright 2004-2010 SRI International
//  Copyright 1997-2009 Torsten Rohlfing
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

#  include <cmtkTypedArray.h>

#  include <string.h>
#  include <stdio.h>

#  include <dcmtk/dcmdata/dcdeftag.h>
#  include <dcmtk/dcmimgle/didocu.h>
#  include <dcmtk/dcmimgle/diutils.h>

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

ScalarImage* 
DICOM::Read 
( const char *path, const Study* study, const int index )
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
  OFCondition status = fileformat->loadFile( path );
  fileformat->transferEnd();

  if ( !status.good() ) 
    {
    StdErr << "Error: cannot read DICOM file " << path << " (" << status.text() << ")\n";
    throw (0);
    }
  
  DcmDataset *dataset = fileformat->getAndRemoveDataset();

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
  const bool haveRescaleIntercept = (0 != document->getValue( DCM_RescaleIntercept, rescaleIntercept ));
  if ( ! haveRescaleIntercept )
    rescaleIntercept = 0;

  const bool haveRescaleSlope = (0 != document->getValue( DCM_RescaleSlope, rescaleSlope ));
  if ( ! haveRescaleSlope )
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

ScalarImage* 
DICOM::Read( const char*, const Study*, const int ) const
{ return NULL; }

void DICOM::Write ( const char*, const ImageInfo&, const int ) {}

//@}

} // namespace cmtk

#endif // #ifdef CMTK_HAVE_DCMTK
