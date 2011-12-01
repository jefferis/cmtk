/*
//
//  Copyright 2004-2011 SRI International
//
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

#include "cmtkDICOM.h"

#include <Base/cmtkTypedArray.h>

#include <System/cmtkConsole.h>

#ifdef CMTK_USE_DCMTK_JPEG
#  include <djdecode.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#include <string.h>
#include <stdio.h>
#include <ctime>

namespace
cmtk
{

/** \addtogroup IO */
//@{

DICOM::DICOM( const char* path )
{
#ifdef CMTK_USE_DCMTK_JPEG
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
    throw Exception( "Could not create DICOM file format object." );
    }
  
  fileformat->transferInit();
  OFCondition status = fileformat->loadFile( path );
  fileformat->transferEnd();

  if ( !status.good() ) 
    {
    throw Exception( "Cannot read DICOM file.." );
//    StdErr << "Error: cannot read DICOM file " << path << " (" << status.text() << ")\n";
    }
  
  this->m_Dataset = fileformat->getAndRemoveDataset();
  
  if ( !this->m_Dataset ) 
    {
    throw Exception( "File format has NULL dataset." );
    }
  
  this->m_Document = std::auto_ptr<DiDocument>( new DiDocument( this->m_Dataset, this->m_Dataset->getOriginalXfer(), CIF_AcrNemaCompatibility ) );
  if ( ! this->m_Document.get() || ! this->m_Document->good() ) 
    {
    throw Exception( "Could not create document representation." );
    } 
}

const FixedVector<3,int>
DICOM::GetDims() const
{
  FixedVector<3,int> dims( FixedVector<3,int>::Init( 0 ) );

  Uint16 tempUint16 = 1;
  if ( this->Document().getValue( DCM_Columns, tempUint16 ) ) 
    {
    dims[0] = static_cast<int>( tempUint16 );
    }

  if ( this->Document().getValue( DCM_Rows, tempUint16 ) ) 
    {
    dims[1] = static_cast<int>( tempUint16 );
    }
  
  // detect and treat multi-frame files
  if ( ! this->Document().getValue( DCM_NumberOfFrames, tempUint16 ) ) 
    {
    // unlike Rows/Columns, NumberofFrames defaults to 1
    tempUint16 = 1;
    }
  dims[2] = tempUint16;

  return dims;
}

const FixedVector<3,double>
DICOM::GetPixelSize() const
{
  FixedVector<3,double> pixelSize( FixedVector<3,double>::Init( 0.0 ) );
  
  // get calibration from image
  const bool hasPixelSpacing = (this->Document().getValue(DCM_PixelSpacing, pixelSize[0], 0) > 0);
  if ( hasPixelSpacing )
    {
    if (this->Document().getValue(DCM_PixelSpacing, pixelSize[1], 1) < 2) 
      {
      throw Exception( "DICOM file does not have two elements in pixel size tag" );
      }
    } 
  else
    throw Exception( "DICOM file does not specify pixel size" );

  // get slice spacing from multi-slice images.
  if ( ! this->Document().getValue( DCM_SpacingBetweenSlices, pixelSize[2] ) )
    {
    pixelSize[2] = 0;
    }

  return pixelSize;
}    
    
const FixedVector<3,double>
DICOM::GetImageOrigin() const
{
  FixedVector<3,double> imageOrigin( FixedVector<3,double>::Init( 0.0 ) );

  const char *image_position_s = NULL;
  if ( ! this->Document().getValue( DCM_ImagePositionPatient, image_position_s ) ) 
    {
    // ImagePositionPatient tag not present, try ImagePosition instead
#ifdef DCM_ImagePosition
    this->Document().getValue( DCM_ImagePosition, image_position_s );
#else
    this->Document().getValue( DCM_ACR_NEMA_ImagePosition, image_position_s );
#endif
    }
  if ( image_position_s ) 
    {
    double xyz[3];
    if ( 3 == sscanf( image_position_s,"%lf%*c%lf%*c%lf", xyz, xyz+1, xyz+2 ) ) 
      {
      imageOrigin = FixedVector<3,double>( xyz );
      }
    }
  
  return imageOrigin;
}

const FixedVector< 2, FixedVector<3,double> > 
DICOM::GetImageOrientation() const
{
  FixedVector< 2, FixedVector<3,double> > orientation;
  
  orientation[0] = FixedVector<3,double>( FixedVector<3,double>::Init( 0.0 ) );
  orientation[1] = FixedVector<3,double>( FixedVector<3,double>::Init( 0.0 ) );

  orientation[0][0] = 1;
  orientation[1][1] = 1;
  
  const char *image_orientation_s = NULL;
#ifdef DCM_ImageOrientation
  if ( ! this->Document().getValue( DCM_ImageOrientation, image_orientation_s ) )
#else
    if ( ! this->Document().getValue( DCM_ACR_NEMA_ImageOrientation, image_orientation_s ) )
#endif
      {
      // ImageOrientation tag not present, try ImageOrientationPatient instead
      this->Document().getValue( DCM_ImageOrientationPatient, image_orientation_s );
      }
  if ( image_orientation_s ) 
    {
    double dx[3], dy[3];
    if ( 6 == sscanf( image_orientation_s, "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", dx, dx+1, dx+2, dy, dy+1, dy+2 ) ) 
      {
      orientation[0] = ( FixedVector<3,double>( dx ) );
      orientation[1] = ( FixedVector<3,double>( dy ) );
      }
    }

  return orientation;
}

TypedArray::SmartPtr
DICOM::GetPixelDataArray( const size_t pixelDataLength )
{
  DcmElement *delem = NULL;

  unsigned short bitsAllocated = 0;
  if ( ( delem = this->m_Document->search( DCM_BitsAllocated ) ) ) 
    {
    delem->getUint16( bitsAllocated );
    } 
  else
    {
    // No "BitsAllocated" tag; use "BitsStored" instead.
    if ( ( delem = this->m_Document->search( DCM_BitsStored ) ) ) 
      {
      delem->getUint16( bitsAllocated );
      }
    }
    
  bool pixelDataSigned = false;
  Uint16 pixelRepresentation = 0;
  if ( this->m_Document->getValue( DCM_PixelRepresentation, pixelRepresentation ) > 0)
    pixelDataSigned = (pixelRepresentation == 1);
    
  double rescaleIntercept, rescaleSlope;
  const bool haveRescaleIntercept = (0 != this->m_Document->getValue( DCM_RescaleIntercept, rescaleIntercept ));
  if ( ! haveRescaleIntercept )
    rescaleIntercept = 0;
    
  const bool haveRescaleSlope = (0 != this->m_Document->getValue( DCM_RescaleSlope, rescaleSlope ));
  if ( ! haveRescaleSlope )
    rescaleSlope = 1;

  pixelDataSigned = pixelDataSigned || (rescaleIntercept < 0);
    
  Uint16 paddingValue = 0;
  const bool paddingFlag = (this->m_Dataset->findAndGetUint16( DCM_PixelPaddingValue, paddingValue )).good();

  TypedArray::SmartPtr pixelDataArray;
    
#ifdef DCM_VariablePixelData
  delem = this->m_Document->search( DCM_VariablePixelData );
#else
  delem = this->m_Document->search( DCM_ACR_NEMA_2C_VariablePixelData );
#endif
  if (!delem)
    delem = this->m_Document->search( DCM_PixelData );
    
  if (delem) 
    {
    if ( (delem->getTag().getEVR() == EVR_OW) || (bitsAllocated > 8) ) 
      {
      Uint16 *pdata = NULL;
      delem->getUint16Array(pdata);
      if ( pixelDataSigned ) 
	{
	const short paddingShort = static_cast<short>( paddingValue );
	pixelDataArray = TypedArray::Create( TYPE_SHORT, pdata, pixelDataLength, paddingFlag, &paddingShort, Memory::ArrayCXX::DeleteWrapper<short> );
	} 
      else
	{
	const unsigned short paddingUShort = static_cast<unsigned short>( paddingValue );
	pixelDataArray = TypedArray::Create( TYPE_USHORT, pdata, pixelDataLength, paddingFlag, &paddingUShort, Memory::ArrayCXX::DeleteWrapper<unsigned short> );
	}
      } 
    else
      {
      Uint8 *pdata = NULL;
      delem->getUint8Array(pdata);
      if ( pixelDataSigned ) 
	{
	const char paddingChar = static_cast<char>( paddingValue );
	pixelDataArray = TypedArray::Create( TYPE_CHAR, pdata, pixelDataLength, paddingFlag, &paddingChar, Memory::ArrayCXX::DeleteWrapper<char> );
	} 
      else 
	{
	const char paddingByte = static_cast<byte>( paddingValue );
	pixelDataArray = TypedArray::Create( TYPE_BYTE, pdata, pixelDataLength, paddingFlag, &paddingByte, Memory::ArrayCXX::DeleteWrapper<byte> );
	}
      }

    delem->detachValueField();
    }

  if ( ! pixelDataArray ) 
    {
    throw( "Could not read pixel data from DICOM file" );
    }
    
  if ( haveRescaleIntercept || haveRescaleSlope ) 
    {
    pixelDataArray->Rescale( rescaleSlope, rescaleIntercept );
    }

  return pixelDataArray;
}


ScalarImage* 
DICOM::Read 
( const char *path )
{
  ScalarImage* image = NULL;

  Self dicom( path );

  FixedVector<3,int> dims = dicom.GetDims();
  FixedVector<3,double> pixelSize = dicom.GetPixelSize();
  ScalarImage::SpaceVectorType imageOrigin = dicom.GetImageOrigin();

  image = new ScalarImage( dims[0], dims[1], dims[2] );
  image->SetPixelSize( pixelSize[0], pixelSize[1] );
  image->SetFrameToFrameSpacing( pixelSize[2] );
  
  TypedArray::SmartPtr pixelDataArray = dicom.GetPixelDataArray( dims[0] * dims[1] * dims[2] );
  image->SetPixelData( pixelDataArray );

  // now some more manual readings...

    // get original table position from image.
  double sliceLocation = 0;
  if ( ! dicom.Document().getValue( DCM_SliceLocation, sliceLocation ) ) 
    {
#ifdef DCM_Location
    dicom.Document().getValue( DCM_Location, sliceLocation );
#else
    dicom.Document().getValue( DCM_ACR_NEMA_Location, sliceLocation );
#endif
    }
  image->SetImageSlicePosition( sliceLocation );
  image->SetImageOrigin( imageOrigin );
  
  // get original image direction from file.
  FixedVector< 2, FixedVector<3,double> > imageOrientation = dicom.GetImageOrientation();
  image->SetImageDirectionX( imageOrientation[0] );
  image->SetImageDirectionY( imageOrientation[1] );

  return image;
}

//@}

} // namespace cmtk
