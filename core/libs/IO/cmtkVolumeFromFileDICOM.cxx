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
#include <Base/cmtkSurfaceNormal.h>

#include <System/cmtkConsole.h>
#include <System/cmtkException.h>

#include <IO/cmtkFileConstHeader.h>
#include <IO/cmtkDICOM.h>

#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmimgle/didocu.h>
#include <dcmtk/dcmimgle/diutils.h>

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

const UniformVolume::SmartPtr
VolumeFromFile::ReadDICOM( const char *path )
{
  DICOM dicom( path );

  DcmElement *delem = NULL;
  Uint16 tempUint16 = 0;

  int dims[3] = { 0, 0, 0 };
  if ( ( delem = dicom.Document().search( DCM_Rows ) ) ) 
    {
    delem->getUint16(tempUint16);
    dims[1]=(int)tempUint16;
    }
    
  if ( ( delem = dicom.Document().search( DCM_Columns ) ) ) 
    {
    delem->getUint16(tempUint16);
    dims[0]=(int)tempUint16;
    }

  // detect and treat multi-frame files
  if ( ! dicom.Document().getValue( DCM_NumberOfFrames, tempUint16 ) ) 
    {
    tempUint16 = 1;
    }
  dims[2] = tempUint16;

  Types::Coordinate pixelSize[3];

  // get calibration from image
  const bool hasPixelSpacing = (dicom.Document().getValue(DCM_PixelSpacing, pixelSize[0], 0) > 0);
  if ( hasPixelSpacing )
    {
    if (dicom.Document().getValue(DCM_PixelSpacing, pixelSize[1], 1) < 2) 
      {
      throw Exception( "DICOM file does not have two elements in pixel size tag" );
      }
    } 
  else
    throw Exception( "DICOM file does not specify pixel size" );
    
  const unsigned long totalImageSizePixels = dims[0] * dims[1] * dims[2];

  TypedArray::SmartPtr pixelDataArray = dicom.GetPixelDataArray();
    

  // now some more manual readings...
    
  // get slice spacing from multi-slice images.
  if ( ! dicom.Document().getValue( DCM_SpacingBetweenSlices, pixelSize[2] ) )
    {
    pixelSize[2] = 0;
    }
    
  // get original image position from file.
  UniformVolume::CoordinateVectorType imageOrigin( UniformVolume::CoordinateVectorType::Init( 0 ) );
  const char *image_position_s = NULL;
  if ( ! dicom.Document().getValue( DCM_ImagePositionPatient, image_position_s ) ) 
    {
    // ImagePositionPatient tag not present, try ImagePosition instead
#ifdef DCM_ImagePosition
    dicom.Document().getValue( DCM_ImagePosition, image_position_s );
#else
    dicom.Document().getValue( DCM_ACR_NEMA_ImagePosition, image_position_s );
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
  if ( ! dicom.Document().getValue( DCM_ImageOrientation, image_orientation_s ) )
#else
    if ( ! dicom.Document().getValue( DCM_ACR_NEMA_ImageOrientation, image_orientation_s ) )
#endif
      {
      // ImageOrientation tag not present, try ImageOrientationPatient
      // instead
      dicom.Document().getValue( DCM_ImageOrientationPatient, image_orientation_s );
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

  // without further information, we "guess" the image normal vector
  UniformVolume::CoordinateVectorType sliceNormal = SurfaceNormal( imageOrientationX, imageOrientationY ).Get();

  // detect and treat Siemens multi-slice mosaics
  const char* tmpStr = NULL;
  if ( dicom.Document().getValue( DCM_Manufacturer, tmpStr ) )
    {
    if ( !strncmp( tmpStr, "SIEMENS", 7 ) )
      {
      const DcmTagKey nSlicesTag(0x0019,0x100a);
      if ( dicom.Document().getValue( nSlicesTag, tempUint16 ) )
	{
	dims[2] = tempUint16;
	
	const DcmTagKey mosaicTag(0x0051,0x100b);
	if ( dicom.Document().getValue( mosaicTag, tmpStr ) )
	  {
	  int rows;
	  int cols;
	  if ( 2 != sscanf( tmpStr, "%dp*%ds", &rows, &cols) )
	    {
	    if ( 2 != sscanf( tmpStr, "%d*%ds", &rows, &cols) )
	      {
	      StdErr << "ERROR: unable to parse mosaic size from " << tmpStr << "\n";
	      }
	    }
	  
	  if ( (cols > 0) && (rows > 0 ) )
	    {
	    const int xMosaic = dims[0] / cols;
	    
	    dims[0] = cols;
	    dims[1] = rows;
	    
	    // de-mosaic the data array
	    const unsigned long imageSizePixels = dims[0] * dims[1] * dims[2];
	    TypedArray::SmartPtr newDataArray( TypedArray::Create( pixelDataArray->GetType(), imageSizePixels ) );
	    
	    const size_t pixelsPerSlice = cols * rows;
	    size_t toOffset = 0;
	    for ( int slice = 0; slice < dims[2]; ++slice )
	      {
	      for ( int j = 0; j < rows; ++j, toOffset += dims[0] )
		{
		const size_t iPatch = slice % xMosaic;
		const size_t jPatch = slice / xMosaic;
		
		const size_t fromOffset = jPatch * xMosaic * pixelsPerSlice + j * xMosaic * cols + iPatch * cols;
		pixelDataArray->BlockCopy( *newDataArray, toOffset, fromOffset, cols );
		}
	      }
	    
	    pixelDataArray = newDataArray;
	    }
	  
// For the following, see here: http://nipy.sourceforge.net/nibabel/dicom/siemens_csa.html#csa-header
	  const DcmTagKey csaHeaderInfoTag(0x0029,0x1010);
	  
	  const Uint8* csaHeaderInfo = NULL;
	  unsigned long csaHeaderLength = 0;
	  dicom.Dataset().findAndGetUint8Array ( csaHeaderInfoTag, csaHeaderInfo, &csaHeaderLength );
	  
	  FileConstHeader fileHeader( csaHeaderInfo, false /*isBigEndian*/ ); // Siemens CSA header is always little endian
	  const size_t nTags = fileHeader.GetField<Uint32>( 8 );
	  
	  size_t tagOffset = 16; // start after header
	  for ( size_t tag = 0; tag < nTags; ++tag )
	    {
	    char tagName[65];
	    fileHeader.GetFieldString( tagOffset, tagName, 64 );
//	  StdErr << tag << "\t" << tagName << "\n";
	    
	    const size_t nItems = fileHeader.GetField<Uint32>( tagOffset + 76 );
//	  StdErr << "  nItems: " << nItems << "\n";
	    
	    tagOffset += 84;
	    for ( size_t item = 0; item < nItems; ++item )
	      {
	      const size_t itemLen = fileHeader.GetField<Uint32>( tagOffset );
	      
//	    StdErr << "    len: " << itemLen << "\n";
	      
	      if ( ! strcmp( tagName, "SliceNormalVector" ) && (item < 3) )
		{
		char valStr[65];
		sliceNormal[item] = atof( fileHeader.GetFieldString( tagOffset+16, valStr, 64 ) );
		
//	      StdErr << "    " << valStr << "\n";
		}
	      
	      tagOffset += 4*((itemLen+3)/4) /*move up to nearest 4-byte boundary*/ + 16 /*the 4 ints at the beginning of item, including itemLength*/;
	      }
	    }
	  }
	}
      }  
    }

  UniformVolume::SmartPtr volume( new UniformVolume( UniformVolume::IndexType( dims ), pixelSize[0], pixelSize[1], pixelSize[2], pixelDataArray ) );
  volume->SetMetaInfo( META_SPACE, "LPS" );
  volume->SetMetaInfo( META_SPACE_ORIGINAL, "LPS" );

  imageOrientationX *= pixelSize[0] / imageOrientationX.RootSumOfSquares();
  imageOrientationY *= pixelSize[1] / imageOrientationY.RootSumOfSquares();
  sliceNormal *= pixelSize[2] / sliceNormal.RootSumOfSquares();

  const Types::Coordinate directions[3][3] = 
    {
      { imageOrientationX[0], imageOrientationX[1], imageOrientationX[2] },
      { imageOrientationY[0], imageOrientationY[1], imageOrientationY[2] },
      { sliceNormal[0], sliceNormal[1], sliceNormal[2] }
    };
  
  const Matrix3x3<Types::Coordinate> m3( directions );
  Matrix4x4<Types::Coordinate> m4( m3 );
  for ( int i = 0; i < 3; ++i )
    m4[3][i] = imageOrigin[i];

  volume->m_IndexToPhysicalMatrix = m4;
  const std::string orientationString0 = volume->GetOrientationFromDirections();
  volume->ChangeCoordinateSpace( AnatomicalOrientation::ORIENTATION_STANDARD );

  const std::string orientationString = volume->GetOrientationFromDirections();
  volume->SetMetaInfo( META_SPACE_UNITS_STRING, "mm" ); // seems to be implied in DICOM
  volume->SetMetaInfo( META_IMAGE_ORIENTATION, orientationString );
  volume->SetMetaInfo( META_IMAGE_ORIENTATION_ORIGINAL, orientationString );

  return volume;
}

} // namespace cmtk
