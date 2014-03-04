/*
//
//  Copyright 1997-2011 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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
#include "nifti1_io_math.h"

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

/// Helper function - make qform in NIFTI header from affine matrix
void
__matrixToNiftiQform( nifti_1_header& header, const AffineXform::MatrixType matrix )
{  
  // transpose and put into nifti_io's mat44 struct
  mat44 R;
  for ( int j = 0; j < 4; ++j )
    {
    for ( int i = 0; i < 4; ++i )
      {
      R.m[i][j] = matrix[j][i];
      }
    }

  // compute quaternion representation
  float qb, qc, qd, qx, qy, qz, dx, dy, dz, qfac;
  nifti_mat44_to_quatern( R, &qb, &qc, &qd, &qx, &qy, &qz, &dx, &dy, &dz, &qfac ) ;
  
  // set header fields
  header.pixdim[0] = qfac;
  header.quatern_b = qb;
  header.quatern_c = qc;
  header.quatern_d = qd;
  header.qoffset_x = qx;
  header.qoffset_y = qy;
  header.qoffset_z = qz;
}

/**\note If the imported NIFTI file contains a qform and/or sform, and the resulting image object is later written back to a new NIFTI file (via VolumeFromFile::WriteNifti), then these transformations will be stored in the output file's qform and sform fields, 
  * respectively.
  * \note If the header contains only a qform, then the qform will be used to initialize cmtk::UniformVolume::m_IndexToPhysicalMatrix, which determines the anatomy-based gross alignment of the image, regardless of qform_code
  * \note If the header contains only an sform, then the sform will be used instead, regardless of sform_code.
  * \note If the header contains both sform and qform, then the transformation is used whose code is NIFTI_XFORM_SCANNER_ANAT. If both sform_code and qform_code are NIFTI_XFORM_SCANNER_ANAT, then qform is used.
  * \note If the header contains neither a qform or an sform, then it is assumed that the image is oriented "RAS", i.e., fastest-moving array index from Left to Right, second fastest from Posterior to Anterior, and third fastest from Inferior to Superior.
  *\warning In keeping with the nifti1.h documentation, it is assumed that qform_code and sform_code will always be non-negative. If this assumption is violated, qform and sform will be switched in the output file (as written by VolumeFromFile::WriteNifti), or
  *   one of them may be missing.
  */
const UniformVolume::SmartPtr
VolumeFromFile::ReadNifti( const std::string& pathHdr, const bool detached, const bool readData )
{
  CompressedStream hdrStream( pathHdr );
  if ( !hdrStream.IsValid() ) 
    {
    StdErr << "ERROR: could not open Nifti header file " << pathHdr << "\n";
    return UniformVolume::SmartPtr( NULL );
    }
  
  nifti_1_header buffer;
  if ( sizeof(buffer) != hdrStream.Read( &buffer, 1, sizeof(buffer) ) ) 
    {
    StdErr << "ERROR: could not read " << sizeof( buffer ) << " bytes from header file " << pathHdr << "\n";
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
    StdErr << "ERROR: image dimension " << ndims << " is smaller than 3 in file " << pathHdr << "\n";
    return UniformVolume::SmartPtr( NULL );
    }
  
  const DataGrid::IndexType::ValueType dims[3] = { header.GetField<short>( 42 ), header.GetField<short>( 44 ), header.GetField<short>( 46 ) };
  const int dims3 = header.GetField<short>( 48 );

  if ( (ndims > 3) && (dims3 > 1) ) 
    {
    StdErr << "WARNING: dimension " << ndims << " is greater than 3 in file " << pathHdr << "\n";
    }
  
  float pixelDim[3];
  header.GetArray( pixelDim, 80, 3 );

  UniformVolume::SmartPtr volume( new UniformVolume( DataGrid::IndexType::FromPointer( dims ), fabs( pixelDim[0] ), fabs( pixelDim[1] ), fabs( pixelDim[2] ) ) );
  // Nifti is in RAS space.
  const char *const niftiSpace = "RAS";
  volume->SetMetaInfo( META_SPACE, niftiSpace );
  volume->SetMetaInfo( META_SPACE_ORIGINAL, niftiSpace );

  const short qform_code = header.GetField<short>( offsetof(nifti_1_header,qform_code) );
  const short sform_code = header.GetField<short>( offsetof(nifti_1_header,sform_code) );

  if ( sform_code > NIFTI_XFORM_UNKNOWN )
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
    if ( (sform_code == NIFTI_XFORM_SCANNER_ANAT) || (qform_code == NIFTI_XFORM_UNKNOWN) )
      volume->m_IndexToPhysicalMatrix = m4;

    // NIFTI says, sform_code must be > 0, so should be safe to use "abs" and make negative to distinguish from qform (below)
    volume->m_AlternativeIndexToPhysicalMatrices[-abs(sform_code)] = m4;
    }
  
  if ( qform_code > NIFTI_XFORM_UNKNOWN )
    {
    const float qb = header.GetField<float>( offsetof(nifti_1_header,quatern_b) );
    const float qc = header.GetField<float>( offsetof(nifti_1_header,quatern_c) );
    const float qd = header.GetField<float>( offsetof(nifti_1_header,quatern_d) );

    const float qx = header.GetField<float>( offsetof(nifti_1_header,qoffset_x) );
    const float qy = header.GetField<float>( offsetof(nifti_1_header,qoffset_y) );
    const float qz = header.GetField<float>( offsetof(nifti_1_header,qoffset_z) );

    const double qfac = (header.GetField<float>( offsetof(nifti_1_header,pixdim) ) >= 0) ? 1.0f : -1.0f;
    const mat44 m44 = nifti_quatern_to_mat44( qb, qc, qd, qx, qy, qz, pixelDim[0], pixelDim[1], pixelDim[2], qfac );

    Matrix4x4<Types::Coordinate> m4;
    for ( int j = 0; j < 4; ++j )
      {
      for ( int i = 0; i < 4; ++i )
	{
	m4[i][j] = m44.m[j][i];
	}
      }

    // qform overrides sform if both exist
    if ( (qform_code == NIFTI_XFORM_SCANNER_ANAT) || (sform_code != NIFTI_XFORM_SCANNER_ANAT) )
      volume->m_IndexToPhysicalMatrix = m4;
    // NIFTI says, qform_code must be > 0, so should be safe to use "abs"
    volume->m_AlternativeIndexToPhysicalMatrices[abs(qform_code)] = m4;
    }

  char orientationImage[4];
  AnatomicalOrientation::GetOrientationFromDirections( orientationImage, volume->m_IndexToPhysicalMatrix, niftiSpace );
  volume->SetMetaInfo( META_IMAGE_ORIENTATION, orientationImage );
  volume->SetMetaInfo( META_IMAGE_ORIENTATION_ORIGINAL, orientationImage );
  
  if ( header.GetField<char>( 148 ) )
    {
    char desc[81];
    desc[80] = 0;
    volume->SetMetaInfo( META_IMAGE_DESCRIPTION, std::string( header.GetFieldString( 148, desc, 80 ) ) );
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
      StdErr << "ERROR: unsupported data type " << header.GetField<short>( 70 ) << " in Nifti file " << pathHdr << "\n";
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

  std::string pathImg = pathHdr;
  if ( detached )
    {
    const size_t period = pathImg.rfind( ".hdr" );
    if ( period != std::string::npos ) 
      pathImg.replace( period, 4, ".img" );
    
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
  
  return volume;
}

void
VolumeFromFile::WriteNifti
( const std::string& path, const UniformVolume& volume )
{
  bool detachedHeader = false;
  bool forceCompressed = false;

  std::string pathImg( path );

  // first, look for .gz
  size_t suffixPosGz = pathImg.rfind( std::string( ".gz" ) );
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

  header.pixdim[0] = 1.0;
  header.pixdim[1] = static_cast<float>( writeVolume->m_Delta[AXIS_X] );
  header.pixdim[2] = static_cast<float>( writeVolume->m_Delta[AXIS_Y] );
  header.pixdim[3] = static_cast<float>( writeVolume->m_Delta[AXIS_Z] );
  header.pixdim[4] = 0.0;
  header.pixdim[5] = 0.0;
  
  for ( std::map<int,cmtk::AffineXform::MatrixType>::const_iterator it = volume.m_AlternativeIndexToPhysicalMatrices.begin(); it != volume.m_AlternativeIndexToPhysicalMatrices.end(); ++it )
    {
    const AffineXform::MatrixType m4 = it->second;

    // this came from a qform
    if ( it->first > 0 ) 
      {
      header.qform_code = it->first;
      __matrixToNiftiQform( header, it->second );
      }
    
    // this came from an sform
    if ( it->first < 0 ) 
      {
      header.sform_code = abs( it->first );
      for ( int i = 0; i < 4; ++i )
	{
	header.srow_x[i] = static_cast<float>( m4[i][0] );
	header.srow_y[i] = static_cast<float>( m4[i][1] );
	header.srow_z[i] = static_cast<float>( m4[i][2] );
	}
      }
    }

  // fallback - we want at least a generic qform to be set to the volume's index-to-physical matrix
  if ( ! (header.qform_code || header.sform_code) )
    {
#ifdef IGNORE
    // This piece of code, when replacing the next two lines (qform stuff) would make CMTK entirely backward-compatible (release 2.3 and earlier), but non-NIFTI-compliant.
    const AffineXform::MatrixType m4 = volume.m_IndexToPhysicalMatrix;
      header.sform_code = NIFTI_XFORM_SCANNER_ANAT;
      for ( int i = 0; i < 4; ++i )
	{
	header.srow_x[i] = static_cast<float>( m4[i][0] );
	header.srow_y[i] = static_cast<float>( m4[i][1] );
	header.srow_z[i] = static_cast<float>( m4[i][2] );
	}
#else
      header.qform_code = NIFTI_XFORM_SCANNER_ANAT;
      __matrixToNiftiQform( header, volume.m_IndexToPhysicalMatrix );
#endif
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
  
  // determine data range;
  const Types::DataItemRange dataRange = data->GetRange();
  header.cal_max = static_cast<float>( dataRange.m_UpperBound );
  header.cal_min = static_cast<float>( dataRange.m_LowerBound );

  if ( volume.MetaKeyExists( META_IMAGE_DESCRIPTION ) )
    {
    memset( header.descrip, 0, sizeof( header.descrip ) );
    strncpy( header.descrip, volume.GetMetaInfo( META_IMAGE_DESCRIPTION ).c_str(), sizeof( header.descrip )-1 );
    }

  if ( volume.MetaKeyExists( META_IMAGE_SLICEORDER ) )
    {
    const std::string sliceOrder = volume.GetMetaInfo( META_IMAGE_SLICEORDER );
    if ( sliceOrder == META_IMAGE_SLICEORDER_SI )
      header.slice_code = NIFTI_SLICE_SEQ_INC;
    else if ( sliceOrder == META_IMAGE_SLICEORDER_SD )
      header.slice_code = NIFTI_SLICE_SEQ_DEC;
    else if ( sliceOrder == META_IMAGE_SLICEORDER_AI )
      header.slice_code = NIFTI_SLICE_ALT_INC;
    else if ( sliceOrder == META_IMAGE_SLICEORDER_AD )
      header.slice_code = NIFTI_SLICE_ALT_DEC;
    else if ( sliceOrder == META_IMAGE_SLICEORDER_AI2 )
      header.slice_code = NIFTI_SLICE_ALT_INC2;
    else if ( sliceOrder == META_IMAGE_SLICEORDER_AD2 )
      header.slice_code = NIFTI_SLICE_ALT_DEC2;

    header.slice_start = 0;
    header.slice_end = header.dim[3]-1; // here we assume that this is a volume straight out of the DICOM stacker, which keeps slices along dim[3]

    header.slice_duration = atof( volume.GetMetaInfo( META_IMAGE_SLICEDURATION ).c_str() );
    }
  
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
      if ( dataSize != CompressedStream::Zlib::StaticSafeWrite( imgFile, data->GetDataPtr(), dataSize ) )
	{
	StdErr << "WARNING: gzwrite() returned error when writing to " << pathImg << "\n";
	}
      gzclose( imgFile );
      }
    else
      {
      StdErr << "ERROR: could not open file '" << pathImg << ".gz' for writing\n";
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
    else
      {
      StdErr << "ERROR: could not open file '" << pathImg << "' for writing\n";
      }
    }
}

} // namespace cmtk
