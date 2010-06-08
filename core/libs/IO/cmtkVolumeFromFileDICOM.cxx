/*
//
//  Copyright 2004-2010 SRI International
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include <cmtkVolumeFromFile.h>


#include <cmtkTypedArray.h>

#include <string.h>
#include <stdio.h>

#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmimgle/didocu.h>
#include <dcmtk/dcmimgle/diutils.h>

#ifdef CMTK_HAVE_DCMTK_JPEG
#  include <djdecode.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_TIME_H
#  include <time.h>
#endif

#include <memory>

#include <cmtkConsole.h>
#include <cmtkException.h>

namespace
cmtk
{

const UniformVolume::SmartPtr
VolumeFromFile::ReadDICOM( const char *path )
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
    throw Exception( "Could not create DICOM file format object." );
    }
    
  fileformat->transferInit();
  fileformat->loadFile( path );
  fileformat->transferEnd();
    
  DcmDataset *dataset = fileformat->getAndRemoveDataset();
  if ( !dataset ) 
    {
    throw Exception( "File format has NULL dataset." );
    }
    
  const E_TransferSyntax xfer = dataset->getOriginalXfer();
  std::auto_ptr<DiDocument> document( new DiDocument( dataset, xfer, CIF_AcrNemaCompatibility ) );
  if ( ! document.get() || ! document->good() ) 
    {
    throw Exception( "Could not create document representation." );
    }
    
  DcmElement *delem = NULL;
  Uint16 tempUint16 = 0;

  int dims[3] = { 0, 0, 0 };
  if ( ( delem = document->search( DCM_Rows ) ) ) 
    {
    delem->getUint16(tempUint16);
    dims[1]=(int)tempUint16;
    }
    
  if ( ( delem = document->search( DCM_Columns ) ) ) 
    {
    delem->getUint16(tempUint16);
    dims[0]=(int)tempUint16;
    }
    
  if ( ! document->getValue( DCM_NumberOfFrames, tempUint16 ) ) 
    {
    tempUint16 = 1;
    }
  dims[2] = tempUint16;

  unsigned short bytesPerPixel = 0;
  unsigned short bitsAllocated = 0;

  if ( ( delem = document->search( DCM_BitsAllocated ) ) ) 
    {
    delem->getUint16( bitsAllocated );
    bytesPerPixel = (((int)bitsAllocated+7)/8);
    } 
  else
    {
    // No "BitsAllocated" tag; use "BitsStored" instead.
    if ( ( delem = document->search( DCM_BitsStored ) ) ) 
      {
      delem->getUint16( bitsAllocated );
      bytesPerPixel = (((int)bitsAllocated+7)/8);
      }
    }
    
  Types::Coordinate pixelSize[3];

  // get calibration from image
  const bool hasPixelSpacing = (document->getValue(DCM_PixelSpacing, pixelSize[0], 0) > 0);
  if ( hasPixelSpacing )
    {
    if (document->getValue(DCM_PixelSpacing, pixelSize[1], 1) < 2) 
      {
      throw Exception( "DICOM file does not have two elements in pixel size tag" );
      }
    } 
  else
    throw Exception( "DICOM file does not specify pixel size" );
    
  const unsigned long totalImageSizePixels = dims[0] * dims[1] * dims[2];
    
  double rescaleIntercept, rescaleSlope;
  const bool haveRescaleIntercept = (0 != document->getValue( DCM_RescaleIntercept, rescaleIntercept ));
  if ( ! haveRescaleIntercept )
    rescaleIntercept = 0;
    
  const bool haveRescaleSlope = (0 != document->getValue( DCM_RescaleSlope, rescaleSlope ));
  if ( ! haveRescaleSlope )
    rescaleSlope = 1;
    
  bool paddingFlag = false;
  Uint16 paddingValue = 0;
  if ( (dataset->findAndGetUint16( DCM_PixelPaddingValue, paddingValue )).good() ) 
    {
    paddingFlag = true;
    } 

  TypedArray::SmartPtr dataArray;
    
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
      dataArray = TypedArray::SmartPtr( TypedArray::Create( TYPE_SHORT, pdata, totalImageSizePixels, true /*freeArray*/ ) );
      } 
    else 
      {
      Uint8 *pdata = NULL;
      delem->getUint8Array(pdata);
      dataArray = TypedArray::SmartPtr( TypedArray::Create( TYPE_CHAR, pdata, totalImageSizePixels, true /*freeArray*/ ) );
      }
      
    if ( paddingFlag )
      dataArray->SetPaddingValue( paddingValue );
      
    if ( haveRescaleIntercept || haveRescaleSlope ) 
      {
      dataArray->Rescale( rescaleSlope, rescaleIntercept );
      }
    delem->detachValueField();
    }
  else
    {
    throw( "Could not read pixel data from DICOM file" );
    }
    
  // now some more manual readings...
    
  // get slice spacing from multi-slice images.
  if ( ! document->getValue( DCM_SpacingBetweenSlices, pixelSize[2] ) )
    {
    pixelSize[2] = 0;
    }
    
  // get original image position from file.
  UniformVolume::CoordinateVectorType imageOrigin( UniformVolume::CoordinateVectorType::Init( 0 ) );
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
    double xyz[3];
    if ( 3 == sscanf( image_position_s,"%lf%*c%lf%*c%lf", xyz, xyz+1, xyz+2 ) ) 
      {
      imageOrigin = UniformVolume::CoordinateVectorType( xyz );
      }
    }
    
  // get original image direction from file.
  UniformVolume::CoordinateVectorType imageOrientationX( UniformVolume::CoordinateVectorType::Init( 0 ) );
  imageOrientationX[0] = 1;
  UniformVolume::CoordinateVectorType imageOrientationY( UniformVolume::CoordinateVectorType::Init( 0 ) );
  imageOrientationY[1] = 1;

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
    double dx[3], dy[3];
    if ( 6 == sscanf( image_orientation_s, "%lf%*c%lf%*c%lf%*c%lf%*c%lf%*c%lf", dx, dx+1, dx+2, dy, dy+1, dy+2 ) ) 
      {
      imageOrientationX = UniformVolume::CoordinateVectorType( dx );
      imageOrientationY = UniformVolume::CoordinateVectorType( dy );
      }
    }
  const Types::Coordinate size[3] = { pixelSize[0] * (dims[0]-1), pixelSize[1] * (dims[1]-1), pixelSize[1] * (dims[1]-1) };  
  return UniformVolume::SmartPtr( new UniformVolume( UniformVolume::IndexType( dims ), UniformVolume::CoordinateVectorType( size ), dataArray ) );
}

} // namespace cmtk
