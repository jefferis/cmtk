/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkScalarImageIO.h>

#include <stdio.h>

#include <cmtkTypes.h>
#include <cmtkFileHeader.h>
#include <cmtkCompressedStream.h>
#include <cmtkAnalyze.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

ScalarImage* 
ScalarImageIO::ReadAnalyze( const char* pathHdr )
{
#ifdef _MSC_VER
  FILE *hdrFile = fopen( pathHdr, "rb" );
#else
  FILE *hdrFile = fopen( pathHdr, "r" );
#endif
  if ( !hdrFile ) 
    {
    StdErr.printf( "ERROR: could not open Analyze header file %s\n",
		      pathHdr );
    return NULL;
    }
  
  char buffer[348];
  if ( 348 != fread( buffer, 1, 348, hdrFile ) ) 
    {
    StdErr.printf( "ERROR: could not read 348 bytes from header file %s\n", pathHdr );
    fclose( hdrFile );
    return NULL;
    }
  fclose( hdrFile );

  const bool bigEndian = (buffer[3] == '\x5C');
  FileHeader header( buffer, bigEndian );

  short ndims = header.GetField<short>( 40 );
  if ( ndims < 2 ) 
    {
    StdErr.printf( "ERROR: image dimension %d is smaller than 2 in file "
		      "%s\n", ndims, pathHdr );
    return NULL;
    }
  
  int dims[4];
  dims[0] = header.GetField<short>( 42 );
  dims[1] = header.GetField<short>( 44 );
  dims[2] = header.GetField<short>( 46 );
  dims[3] = header.GetField<short>( 48 );

  if ( (ndims > 2) && ((dims[2] > 1) || (dims[3] > 1)) ) 
    {
    StdErr.printf( "WARNING: dimension %d is greater than 2 in file %s\n",
		      ndims, pathHdr );
    }
  
  float pixelDim[2];
  header.GetArray( pixelDim, 80, 2 );

  ScalarImage* image = new ScalarImage( dims[0], dims[1] );
  image->SetPixelSize( pixelDim[0], pixelDim[1] );

  ScalarDataType dtype;
  switch ( header.GetField<short>( 70 ) ) 
    {
    case cmtk::ANALYZE_TYPE_NONE:
    case cmtk::ANALYZE_TYPE_BINARY:
    case cmtk::ANALYZE_TYPE_COMPLEX:
    case cmtk::ANALYZE_TYPE_RGB:
    case cmtk::ANALYZE_TYPE_ALL:
    default:
    StdErr.printf( "ERROR: unsupported data type in Analyze file %s\n",
		      pathHdr );
    return NULL;
    case cmtk::ANALYZE_TYPE_UNSIGNED_CHAR:
      dtype = TYPE_BYTE;
      break;
    case cmtk::ANALYZE_TYPE_SIGNED_SHORT:
      dtype = TYPE_SHORT;
      break;
    case cmtk::ANALYZE_TYPE_SIGNED_INT:
      dtype = TYPE_INT;
      break;
    case cmtk::ANALYZE_TYPE_FLOAT:
      dtype = TYPE_FLOAT;
      break;
    case cmtk::ANALYZE_TYPE_DOUBLE:
      dtype = TYPE_DOUBLE;
      break;
    }  
  image->CreatePixelData( dtype );
  
  size_t offset = static_cast<size_t>( header.GetField<float>( 108 ) );

  char* pathImg = Memory::AllocateArray<char>(  4 + strlen( pathHdr )  );
  strcpy( pathImg, pathHdr );
  char* suffix = strstr( pathImg, ".hdr" );
  if ( suffix ) *suffix = 0;
  strcat( pathImg, ".img" );

  CompressedStream stream( pathImg );
  if ( stream.IsValid() ) 
    {
    stream.Seek( offset, SEEK_CUR );
    
    TypedArray::SmartPtr data = image->GetPixelData();
    stream.Read
      ( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize() );
    
#ifdef WORDS_BIGENDIAN
    if ( ! bigEndian ) data->ChangeEndianness();
#else
    if ( bigEndian ) data->ChangeEndianness();
#endif
    
    } 
  else 
    {
    StdErr.printf( "WARNING: could not open Analyze image file %s\n",
		      pathImg );
    }
  
  delete[] pathImg;

  return image;
}

void
ScalarImageIO::WriteAnalyze
( const char* pathHdr, const ScalarImage* image )
{
  const TypedArray* data = image->GetPixelData();
  if ( ! data ) return;

  char buffer[348];
#ifdef WORDS_BIGENDIAN
  FileHeader header( buffer, true /*bigEndian*/ );
#else
  FileHeader header( buffer, false /*bigEndian*/ );
#endif

  header.StoreField<int>( 0, 348 ); // header size
  header.StoreField<int>( 32, 16384 ); // extents
  header.StoreField<short>( 36, 0 ); // session error
  header.StoreField<char>( 38, 'r' ); // regular

  // ndims
  header.StoreField<short>( 40, 4 );

  // dimensions
  header.StoreField<short>( 42, image->GetDims( AXIS_X ) );
  header.StoreField<short>( 44, image->GetDims( AXIS_Y ) );
  header.StoreField<short>( 46, 1 );
  header.StoreField<short>( 48, 1 ); // write dims 3-7
  header.StoreField<short>( 50, 0 ); // just for safety
  header.StoreField<short>( 52, 0 ); // just for safety
  header.StoreField<short>( 54, 0 ); // just for safety

  header.StoreField<float>( 68, 0.0 ); // vox_offset
  switch ( data->GetType() ) 
    {
    default:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_NONE );
      header.StoreField<short>( 72, 0 );
      break;
    case TYPE_BYTE:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_UNSIGNED_CHAR );
      header.StoreField<short>( 72, 8 * sizeof( unsigned char ) );
      break;
    case TYPE_SHORT:
    case TYPE_USHORT:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_SIGNED_SHORT );
      header.StoreField<short>( 72, 8 * sizeof( short ) );
      break;
    case TYPE_INT:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_SIGNED_INT );
      header.StoreField<short>( 72, 8 * sizeof( int ) );
      break;
    case TYPE_FLOAT:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_FLOAT );
      header.StoreField<short>( 72, 8 * sizeof( float ) );
      break;
    case TYPE_DOUBLE:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_DOUBLE );
      header.StoreField<short>( 72, 8 * sizeof( double ) );
      break;
    }  
  
  header.StoreField<float>( 80, (float)image->GetPixelSize( AXIS_X ) );
  header.StoreField<float>( 84, (float)image->GetPixelSize( AXIS_Y ) );
  header.StoreField<float>( 88, 1.0f );
  header.StoreField<float>( 92, 1.0f ); // write sizes in dims 3 and
  header.StoreField<float>( 96, 1.0f ); // 4 just to be safe

  // set zero offset for binary file.
  header.StoreField<float>( 108, 0.0f );

  // determine data range;
  Types::DataItem dataMin, dataMax;
  data->GetRange( dataMin, dataMax );

  header.StoreField<float>( 124, static_cast<float>( dataMax ) ); // cal_max
  header.StoreField<float>( 128, static_cast<float>( dataMin ) ); // cal_min

  header.StoreField<int>( 140, static_cast<int>( dataMax ) );
  header.StoreField<int>( 144, static_cast<int>( dataMin ) );

  // slice orientation always axial from caudal to cranial
  header.StoreField<byte>( 252, 0 );

  // write header info
#ifdef MSC_VER
  FILE *hdrFile = fopen( pathHdr, "wb" );
#else
  FILE *hdrFile = fopen( pathHdr, "w" );
#endif
  if ( hdrFile )
    {
    if ( 348 != fwrite( buffer, 1, 348, hdrFile ) ) 
      {
      StdErr.printf( "ERROR: could not write 348 bytes to header file %s\n",
			pathHdr );
      }
    fclose( hdrFile );
    }

  // write binary data
  char* pathImg = Memory::AllocateArray<char>(  4 + strlen( pathHdr )  );
  strcpy( pathImg, pathHdr );
  char* suffix = strstr( pathImg, ".hdr" );
  if ( suffix ) *suffix = 0;
  strcat( pathImg, ".img" );

#ifdef _MSC_VER
  FILE *imgFile = fopen( pathImg, "wb" );
#else
  FILE *imgFile = fopen( pathImg, "w" );
#endif
  if ( imgFile ) 
    {
    fwrite( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize(), imgFile );
    fclose( imgFile );
    }

  delete[] pathImg;
}

} // namespace cmtk
