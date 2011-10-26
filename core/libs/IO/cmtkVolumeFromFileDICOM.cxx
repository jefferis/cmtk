/*
//
//  Copyright 2004-2011 SRI International
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

#include "cmtkVolumeFromFile.h"

#include <Base/cmtkTypedArray.h>

#include <System/cmtkConsole.h>
#include <System/cmtkException.h>

#include <IO/cmtkFileConstHeader.h>

#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmimgle/didocu.h>
#include <dcmtk/dcmimgle/diutils.h>

#ifdef CMTK_USE_DCMTK_JPEG
#  include <djdecode.h>
#endif

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#include <memory>
#include <string.h>
#include <stdio.h>
#include <ctime>

namespace
cmtk
{

const UniformVolume::SmartPtr
VolumeFromFile::ReadDICOM( const char *path )
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

  // detect and treat multi-frame files
  if ( ! document->getValue( DCM_NumberOfFrames, tempUint16 ) ) 
    {
    tempUint16 = 1;
    }
  dims[2] = tempUint16;

  unsigned short bitsAllocated = 0;
  if ( ( delem = document->search( DCM_BitsAllocated ) ) ) 
    {
    delem->getUint16( bitsAllocated );
    } 
  else
    {
    // No "BitsAllocated" tag; use "BitsStored" instead.
    if ( ( delem = document->search( DCM_BitsStored ) ) ) 
      {
      delem->getUint16( bitsAllocated );
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

  // detect and treat Siemens multi-slice mosaics
  const char* tmpStr = NULL;
  if ( document->getValue( DCM_Manufacturer, tmpStr ) )
    {
    if ( !strncmp( tmpStr, "SIEMENS", 7 ) )
      {
      const DcmTagKey mosaicTag(0x0051,0x100b);
      
      if ( document->getValue( mosaicTag, tmpStr ) )
	{
	int rows;
	int cols;
	sscanf( tmpStr, "%dp*%ds", &rows, &cols);

	const int xMosaic = dims[0] / cols;
	const int yMosaic = dims[1] / rows;
	
	dims[0] = cols;
	dims[1] = rows;

	const DcmTagKey nSlicesTag(0x0019,0x100a);
	if ( document->getValue( nSlicesTag, tempUint16 ) )
	  {
	  dims[2] = tempUint16;

	  // de-mosaic the data array
	  const unsigned long imageSizePixels = dims[0] * dims[1] * dims[2];
	  TypedArray::SmartPtr newDataArray( TypedArray::Create( dataArray->GetType(), imageSizePixels ) );

	  const size_t pixelsPerSlice = cols * rows;
	  size_t toOffset = 0;
	  for ( int slice = 0; slice < dims[2]; ++slice )
	    {
	    for ( int j = 0; j < rows; ++j, toOffset += dims[0] )
	      {
	      const size_t iPatch = slice % xMosaic;
	      const size_t jPatch = slice / xMosaic;

	      const size_t fromOffset = jPatch * xMosaic * pixelsPerSlice + j * xMosaic * cols + iPatch * cols;
	      dataArray->BlockCopy( *newDataArray, toOffset, fromOffset, cols );
	      }
	    }

	  dataArray = newDataArray;
	  }
	else
	  {
	  StdErr << "WARNING: image claims to by a Siemens mosaic file, but does not provide number of slices in tag (0x0019,0x100a)\n";
	  }
	
// For the following, see here: http://nipy.sourceforge.net/nibabel/dicom/siemens_csa.html#csa-header
	const DcmTagKey csaHeaderInfoTag(0x0029,0x1010);

	const Uint8* csaHeaderInfo = NULL;
	unsigned long csaHeaderLength = 0;
	dataset->findAndGetUint8Array ( csaHeaderInfoTag, csaHeaderInfo, &csaHeaderLength );

	FileConstHeader fileHeader( csaHeaderInfo, false /*isBigEndian*/ ); // Siemens CSA header is always little endian
	const size_t nTags = fileHeader.GetField<Uint32>( 8 );

	size_t tagOffset = 16; // start after header
	for ( int tag = 0; tag < nTags; ++tag )
	  {
	  char tagName[65];
	  fileHeader.GetFieldString( tagOffset, tagName, 64 );
	  StdErr << tag << "\t" << tagName << "\n";

	  const size_t nItems = fileHeader.GetField<Uint32>( tagOffset + 76 );
	  StdErr << "  nItems: " << nItems << "\n";

	  tagOffset += 84;
	  for ( int item = 0; item < nItems; ++item )
	    {
	    const size_t itemLen = fileHeader.GetField<Uint32>( tagOffset );

	    StdErr << "    len: " << itemLen << "\n";
	    tagOffset += 4*((itemLen+3)/4) + 16 /*the 4 ints at the beginning of item, including itemLength*/;
	    }
	  }
	}
      }
    }  

  UniformVolume::SmartPtr volume( new UniformVolume( UniformVolume::IndexType( dims ), pixelSize[0], pixelSize[1], pixelSize[2], dataArray ) );
  volume->SetMetaInfo( META_SPACE, "LPS" );
  volume->SetMetaInfo( META_SPACE_ORIGINAL, "LPS" );

  return volume;
}

} // namespace cmtk
