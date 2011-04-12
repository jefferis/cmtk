/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include <System/cmtkConsole.h>
#include <System/cmtkCompressedStream.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkAnalyze.h>
#include <IO/cmtkFileHeader.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkAnatomicalOrientation.h>

#include "nifti1.h"

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#ifdef HAVE_ZLIB
#  include <zlib.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

const UniformVolume::SmartPtr
VolumeFromFile::ReadNifti( const char* pathHdr, const bool detached, const bool readData )
{
  CompressedStream hdrStream( pathHdr );
  if ( !hdrStream.IsValid() ) 
    {
    StdErr.printf( "ERROR: could not open Nifti header file %s\n", pathHdr );
    return UniformVolume::SmartPtr( NULL );
    }
  
  nifti_1_header buffer;
  if ( sizeof(buffer) != hdrStream.Read( &buffer, 1, sizeof(buffer) ) ) 
    {
    StdErr.printf( "ERROR: could not read %d bytes from header file %s\n", int( sizeof( buffer ) ), pathHdr );
    return UniformVolume::SmartPtr( NULL );
    }
  hdrStream.Close();
  
  // determine if we need to byte-swap
  const int dim0 = buffer.dim[0];
  const bool byteSwap = ((dim0>0) && (dim0<8)) ? false : true;

#ifdef WORDS_BIGENDIAN
  FileHeader header( &buffer, !byteSwap );
#else
  FileHeader header( &buffer, byteSwap );
#endif

  short ndims = header.GetField<short>( 40 );
  if ( ndims < 3 ) 
    {
    StdErr.printf( "ERROR: image dimension %d is smaller than 3 in file %s\n", ndims, pathHdr );
    return UniformVolume::SmartPtr( NULL );
    }
  
  const DataGrid::IndexType::ValueType dims[3] = { header.GetField<short>( 42 ), header.GetField<short>( 44 ), header.GetField<short>( 46 ) };
  const int dims3 = header.GetField<short>( 48 );

  if ( (ndims > 3) && (dims3 > 1) ) 
    {
    StdErr.printf( "WARNING: dimension %d is greater than 3 in file %s\n", ndims, pathHdr );
    }
  
  float pixelDim[3];
  header.GetArray( pixelDim, 80, 3 );

  Types::Coordinate size[3] = { (dims[0] - 1) * fabs( pixelDim[0] ), (dims[1] - 1) * fabs( pixelDim[1] ), (dims[2] - 1) * fabs( pixelDim[2] ) };

  UniformVolume::SmartPtr volume( new UniformVolume( DataGrid::IndexType( dims ), UniformVolume::CoordinateVectorType( size ) ) );
  // Nifti is in RAS space.
  const char *const niftiSpace = "RAS";
  volume->m_MetaInformation[META_SPACE] = volume->m_MetaInformation[META_SPACE_ORIGINAL] = niftiSpace;

  const short qform_code = header.GetField<short>( offsetof(nifti_1_header,qform_code) );
  if ( qform_code > 0 )
    {
    const float qb = header.GetField<float>( offsetof(nifti_1_header,quatern_b) );
    const float qc = header.GetField<float>( offsetof(nifti_1_header,quatern_c) );
    const float qd = header.GetField<float>( offsetof(nifti_1_header,quatern_d) );
    const double qa = sqrt( std::max( 0.0, 1.0 - (qb*qb + qc*qc + qd*qd) ) );

    const float qfac = (header.GetField<float>( offsetof(nifti_1_header,pixdim) ) >= 0) ? 1.0f : -1.0f;

    const Types::Coordinate directions[3][3] = 
      {
	{        pixelDim[0] * (qa*qa+qb*qb-qc*qc-qd*qd),        pixelDim[0] * (2*qb*qc+2*qa*qd),                pixelDim[0] * (2*qb*qd-2*qa*qc) },
	{        pixelDim[1] * (2*qb*qc-2*qa*qd),                pixelDim[1] * (qa*qa+qc*qc-qb*qb-qd*qd),        pixelDim[1] * (2*qc*qd+2*qa*qb) },
	{ qfac * pixelDim[2] * (2*qb*qd+2*qa*qc),         qfac * pixelDim[2] * (2*qc*qd-2*qa*qb),         qfac * pixelDim[2] * (qa*qa+qd*qd-qc*qc-qb*qb) }
      };
    
    const Matrix3x3<Types::Coordinate> m3( directions );
    Matrix4x4<Types::Coordinate> m4( m3 );

    m4[3][0] = header.GetField<float>( offsetof(nifti_1_header,qoffset_x) );
    m4[3][1] = header.GetField<float>( offsetof(nifti_1_header,qoffset_y) );
    m4[3][2] = header.GetField<float>( offsetof(nifti_1_header,qoffset_z) );

    volume->m_IndexToPhysicalMatrix = m4;
    
    char orientationImage[4];
    AnatomicalOrientation::GetOrientationFromDirections( orientationImage, m4, niftiSpace );
    volume->m_MetaInformation[META_IMAGE_ORIENTATION] = volume->m_MetaInformation[META_IMAGE_ORIENTATION_ORIGINAL] = orientationImage;
    }
  else
    {
    const short sform_code = header.GetField<short>( offsetof(nifti_1_header,sform_code) );
    if ( sform_code > 0 )
      {
      float srow_x[4], srow_y[4], srow_z[4];
      header.GetArray( srow_x, offsetof(nifti_1_header,srow_x), 4 );
      header.GetArray( srow_y, offsetof(nifti_1_header,srow_y), 4 );
      header.GetArray( srow_z, offsetof(nifti_1_header,srow_z), 4 );

      const Types::Coordinate directions[4][4] = 	
      {
	{ srow_x[0], srow_y[0], srow_z[0], 0 },
	{ srow_x[1], srow_y[1], srow_z[1], 0 },
	{ srow_x[2], srow_y[2], srow_z[2], 0 },
	{ srow_x[3], srow_y[3], srow_z[3], 1 }
      };
      
      Matrix4x4<Types::Coordinate> m4( directions );      
      volume->m_IndexToPhysicalMatrix = m4;
    
      char orientationImage[4];
      AnatomicalOrientation::GetOrientationFromDirections( orientationImage, m4, niftiSpace );
      volume->m_MetaInformation[META_IMAGE_ORIENTATION] = volume->m_MetaInformation[META_IMAGE_ORIENTATION_ORIGINAL] = orientationImage;
      }
    else
      {
      // no orientation info, default to nifti standard space
      volume->m_MetaInformation[META_IMAGE_ORIENTATION] = volume->m_MetaInformation[META_IMAGE_ORIENTATION_ORIGINAL] = niftiSpace;
      }
    }
  
  // don't read data, we're done here.
  if ( ! readData )
    return volume;

  ScalarDataType dtype;
  switch ( header.GetField<short>( 70 ) ) 
    {
    case DT_NONE:
    case DT_BINARY:
    case DT_COMPLEX:
    case DT_COMPLEX128:
    case DT_COMPLEX256:
    case DT_INT64:
    case DT_UINT64:
    case DT_FLOAT128:
    case DT_RGB:
    case DT_ALL:
    default:
      StdErr.printf( "ERROR: unsupported data type %d in Nifti file %s\n", header.GetField<short>( 70 ), pathHdr );
      return volume;
    case DT_UNSIGNED_CHAR:
      dtype = TYPE_BYTE;
      break;
    case DT_INT8:
      dtype = TYPE_CHAR;
      break;
    case DT_SIGNED_SHORT:
      dtype = TYPE_SHORT;
      break;
    case DT_SIGNED_INT:
      dtype = TYPE_INT;
      break;
    case DT_FLOAT:
      dtype = TYPE_FLOAT;
      break;
    case DT_DOUBLE:
      dtype = TYPE_DOUBLE;
      break;
    case DT_UINT16:
      dtype = TYPE_USHORT;
      break;
    case DT_UINT32:
      dtype = TYPE_UINT;
      break;
    }  
  
  size_t offset = static_cast<size_t>( header.GetField<float>( 108 ) );
  char* pathImg = Memory::AllocateArray<char>(  4 + strlen( pathHdr )  );
  strcpy( pathImg, pathHdr );

  if ( detached )
    {
    char* suffix = strstr( pathImg, ".hdr" );
    if ( suffix ) *suffix = 0;
    strcat( pathImg, ".img" );
    offset = 0;
    }
  
  CompressedStream stream( pathImg );
  if ( stream.IsValid() ) 
    {
    stream.Seek( offset, SEEK_CUR );
    
    TypedArray::SmartPtr data( TypedArray::Create( dtype, volume->GetNumberOfPixels() ) );
    if ( data->GetDataSize() == stream.Read( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize() ) )
      {
      if ( byteSwap ) 
	data->ChangeEndianness();
      
      volume->SetData( data );
      }
    else
      {
      StdErr << "ERROR: could not read " << data->GetDataSize() << " pixels from Nifti image file " << pathImg << "\n";
      }
    } 
  else
    {
    StdErr << "ERROR: could not open Nifti image file " << pathImg << "\n";
    }
  
  Memory::DeleteArray( pathImg );
  
  return volume;
}

void
VolumeFromFile::WriteNifti
( const char* path, const UniformVolume& volume )
{
  bool detachedHeader = false;
  bool forceCompressed = false;

  std::string pathImg( path );

  // first, look for .gz
  size_t suffixPosGz = pathImg.rfind( ".gz" );
  if ( suffixPosGz != std::string::npos )
    {
    // found: set force compression flag and remove .gz from path
    forceCompressed = true;
    pathImg = pathImg.substr( 0, suffixPosGz );
    }
  
  std::string pathHdr( pathImg );
  size_t suffixPos = pathHdr.rfind( ".img" );
  if ( suffixPos != std::string::npos )
    {
    detachedHeader = true;
    pathHdr.replace( suffixPos, 4, ".hdr" );
    }
  
  UniformVolume::SmartPtr writeVolume( volume.Clone() );
  writeVolume->ChangeCoordinateSpace( "RAS" );

  const TypedArray* data = writeVolume->GetData();
  if ( ! data ) return;

  nifti_1_header header;
  memset( &header, 0, sizeof( header ) );

  header.sizeof_hdr = 348; // header size
  header.dim_info = 0;

  // ndims
  header.dim[0] = 4;

  // dimensions
  header.dim[1] = writeVolume->GetDims()[AXIS_X];
  header.dim[2] = writeVolume->GetDims()[AXIS_Y];
  header.dim[3] = writeVolume->GetDims()[AXIS_Z];
  header.dim[4] = 1;
  header.dim[5] = 0;
  header.dim[6] = 0;
  header.dim[7] = 0;

  header.qform_code = 0;
  header.sform_code = 1;

  AffineXform::MatrixType m4 = volume.m_IndexToPhysicalMatrix;
  for ( int i = 0; i < 4; ++i )
    {
    header.srow_x[i] = static_cast<float>( m4[i][0] );
    header.srow_y[i] = static_cast<float>( m4[i][1] );
    header.srow_z[i] = static_cast<float>( m4[i][2] );
    }
  
  switch ( data->GetType() ) 
    {
    default:
      header.datatype = DT_UNKNOWN;
      header.bitpix = 0;
      break;
    case TYPE_BYTE:
      header.datatype = DT_UNSIGNED_CHAR;
      header.bitpix = 8 * sizeof(byte);
      break;
    case TYPE_CHAR:
      header.datatype = DT_INT8;
      header.bitpix = 8 * sizeof(char);
      break;
    case TYPE_SHORT:
      header.datatype = DT_INT16;
      header.bitpix = 8 * sizeof(short);
      break;
    case TYPE_USHORT:
      header.datatype = DT_UINT16;
      header.bitpix = 8 * sizeof(unsigned short);
      break;
    case TYPE_INT:
      header.datatype = DT_INT32;
      header.bitpix = 8 * sizeof(int);
      break;
    case TYPE_UINT:
      header.datatype = DT_UINT32;
      header.bitpix = 8 * sizeof(unsigned int);
      break;
    case TYPE_FLOAT:
      header.datatype = DT_FLOAT;
      header.bitpix = 8 * sizeof(float);
      break;
    case TYPE_DOUBLE:
      header.datatype = DT_DOUBLE;
      header.bitpix = 8 * sizeof(double);
      break;
    }  
  
  header.pixdim[0] = 1;
  header.pixdim[1] = static_cast<float>( writeVolume->m_Delta[AXIS_X] );
  header.pixdim[2] = static_cast<float>( writeVolume->m_Delta[AXIS_Y] );
  header.pixdim[3] = static_cast<float>( writeVolume->m_Delta[AXIS_Z] );
  header.pixdim[4] = 0.0;
  header.pixdim[5] = 0.0;
  
  // determine data range;
  const Types::DataItemRange dataRange = data->GetRange();
  header.cal_max = static_cast<float>( dataRange.m_UpperBound );
  header.cal_min = static_cast<float>( dataRange.m_LowerBound );
  
#ifdef _MSC_VER
  const char *const modestr = "wb";
#else
  const char *const modestr = "w";
#endif
    
  if ( detachedHeader )
    {
    memcpy( &header.magic, "ni1\x00", 4 );
    header.vox_offset = 0;
    FILE *hdrFile = fopen( pathHdr.c_str(), modestr );
    if ( hdrFile ) 
      {
      fwrite( &header, 1, sizeof( header ), hdrFile );
      const int extension = 0;
      fwrite( &extension, 1, 4, hdrFile );
      fclose( hdrFile );
      }
    else
      {
      StdErr << "ERROR: NIFTI header file '" << pathHdr << "' could not be opened for writing!\n";
      }
    }
  else
    {
    memcpy( &header.magic, "n+1\x00", 4 );
    header.vox_offset = 352;
    }
  
  if ( VolumeIO::GetWriteCompressed() || forceCompressed )
    {
    struct stat buf;
    if ( ! stat( pathImg.c_str(), &buf ) )
      {
      StdErr << "WARNING: NIFTI file '" << path << "' will be written compressed, but uncompressed file exists!\n";
      }
    
    gzFile imgFile = gzopen( (pathImg+".gz").c_str(), modestr );
    if ( imgFile ) 
      {
      if ( ! detachedHeader )
	{
	gzwrite( imgFile, &header, sizeof( header ) );
	const int extension = 0;
	gzwrite( imgFile, &extension, 4 );
	}
      
      const size_t dataSize = data->GetItemSize() * data->GetDataSize();
      if ( dataSize != static_cast<size_t>( gzwrite( imgFile, data->GetDataPtr(), dataSize ) ) )
	{
	StdErr << "WARNING: gzwrite() returned error when writing to " << pathImg << "\n";
	}
      gzclose( imgFile );
      }
    }
  else
    {
    FILE *imgFile = fopen( pathImg.c_str(), modestr );
    if ( imgFile ) 
      {
      if ( ! detachedHeader )
	{
	fwrite( &header, 1, sizeof( header ), imgFile );
	const int extension = 0;
	fwrite( &extension, 1, 4, imgFile );
	}

      fwrite( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize(), imgFile );
      fclose( imgFile );
      }
    }
}

} // namespace cmtk
