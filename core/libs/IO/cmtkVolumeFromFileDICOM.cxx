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

  FixedVector<3,int> dims = dicom.GetDims();
  FixedVector<3,double> pixelSize = dicom.GetPixelSize();

  const unsigned long totalImageSizePixels = dims[0] * dims[1] * dims[2];

  TypedArray::SmartPtr pixelDataArray = dicom.GetPixelDataArray( totalImageSizePixels );

  const UniformVolume::CoordinateVectorType imageOrigin = dicom.GetImageOrigin();
  FixedVector< 2, FixedVector<3,double> > imageOrientation = dicom.GetImageOrientation();
  
  // now some more manual readings...
    
  // without further information, we "guess" the image normal vector
  UniformVolume::CoordinateVectorType sliceNormal = SurfaceNormal( imageOrientation[0], imageOrientation[1] ).Get();

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

  imageOrientation[0] *= pixelSize[0] / imageOrientation[0].RootSumOfSquares();
  imageOrientation[1] *= pixelSize[1] / imageOrientation[1].RootSumOfSquares();
  sliceNormal *= pixelSize[2] / sliceNormal.RootSumOfSquares();

  const Types::Coordinate directions[3][3] = 
    {
      { imageOrientation[0][0], imageOrientation[0][1], imageOrientation[0][2] },
      { imageOrientation[1][0], imageOrientation[1][1], imageOrientation[1][2] },
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
