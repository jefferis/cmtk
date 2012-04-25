/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
#include <System/cmtkDebugOutput.h>
#include <System/cmtkCompressedStream.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkAnalyze.h>
#include <IO/cmtkFileHeader.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkAnatomicalOrientation.h>

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

/// Environment variable that turns on legacy Analyze I/O with incorrect anatomical orientations.
const char* const CMTK_LEGACY_ANALYZE_IO = "CMTK_LEGACY_ANALYZE_IO";

/// Legacy environment variable that turns on legacy Analyze I/O with incorrect anatomical orientations.
const char* const IGS_LEGACY_ANALYZE_IO = "IGS_LEGACY_ANALYZE_IO";

/** \addtogroup IO */
//@{

const UniformVolume::SmartPtr
VolumeFromFile::ReadAnalyzeHdr( const std::string& pathHdr, const bool bigEndian, const bool readData )
{
#ifdef _MSC_VER
  FILE *hdrFile = fopen( pathHdr.c_str(), "rb" );
#else
  FILE *hdrFile = fopen( pathHdr.c_str(), "r" );
#endif
  if ( !hdrFile ) 
    {
    StdErr << "ERROR: could not open Analyze header file " << pathHdr << "\n";
    return UniformVolume::SmartPtr( NULL );
    }
  
  char buffer[348];
  if ( 348 != fread( buffer, 1, 348, hdrFile ) ) 
    {
    StdErr << "ERROR: could not read 348 bytes from header file " << pathHdr << "\n";
    fclose( hdrFile );
    return UniformVolume::SmartPtr( NULL );
    }
  fclose( hdrFile );

  FileHeader header( buffer, bigEndian );

  short ndims = header.GetField<short>( 40 );
  if ( ndims < 3 ) 
    {
    StdErr.printf( "ERROR: image dimension %d is smaller than 3 in file %s\n", ndims, pathHdr.c_str() );
    return UniformVolume::SmartPtr( NULL );
    }
  
  DataGrid::IndexType dims;
  dims[0] = header.GetField<short>( 42 );
  dims[1] = header.GetField<short>( 44 );
  dims[2] = header.GetField<short>( 46 );
  const int dims3 = header.GetField<short>( 48 );

  if ( (ndims > 3) && (dims3 > 1) ) 
    {
    StdErr.printf( "WARNING: dimension %d is greater than 3 in file %s\n", ndims, pathHdr.c_str() );
    }
  
  float pixelDim[3];
  header.GetArray( pixelDim, 80, 3 );

// use fabs to catch something weird in FSL's output
#ifndef CMTK_REGRESSION
  const Types::Coordinate size[3] = { float((dims[0] - 1) * fabs( pixelDim[0] )), float((dims[1] - 1) * fabs( pixelDim[1] )), float((dims[2] - 1) * fabs( pixelDim[2] )) };
#else
  const Types::Coordinate size[3] = { (dims[0] - 1) * fabs( pixelDim[0] ), (dims[1] - 1) * fabs( pixelDim[1] ), (dims[2] - 1) * fabs( pixelDim[2] ) };
#endif
  
  const byte orient = header.GetField<byte>( 252 );

  const char* orientString = NULL;

  const bool legacyMode = !header.CompareFieldStringN( 344, "SRI1", 4 );
  if ( legacyMode )
    {
    // 
    // Legacy mode: read Analyze images with incorrect orientations to preserve backward compatibility
    // We do this by setting "apparent" anatomical axis orientation codes so that the new reorientation
    // function in DataGrid will perform the same reordering as its obsolete SetData() member used to.
    //
    switch ( orient ) 
      {
      default:
	StdErr.printf( "WARNING: unsupported slice orientation %d in Analyze file %s\n", orient, pathHdr.c_str() );
      case cmtk::ANALYZE_AXIAL:
	orientString = "RAS"; // INCORRECT LEGACY ORIENTATION
	break;
      case cmtk::ANALYZE_AXIAL_FLIP:
	orientString = "RAI"; // INCORRECT LEGACY ORIENTATION
	break;
      case cmtk::ANALYZE_CORONAL:
	orientString = "RIP"; // INCORRECT LEGACY ORIENTATION
	break;
      case cmtk::ANALYZE_CORONAL_FLIP:
	orientString = "RSA"; // INCORRECT LEGACY ORIENTATION
	break;
      case cmtk::ANALYZE_SAGITTAL:
	orientString = "AIR"; // INCORRECT LEGACY ORIENTATION
	break;
      case cmtk::ANALYZE_SAGITTAL_FLIP:
	orientString = "AIL"; // INCORRECT LEGACY ORIENTATION
	break;
      }
    StdErr << "INFO: reading Analyze hdr/img in legacy orientation mode, assuming " << orientString << " axes\n";
    }
  else
    {
    switch ( orient ) 
      {
      default:
	StdErr.printf( "WARNING: unsupported slice orientation %d in Analyze file %s\n", orient, pathHdr.c_str() );
      case cmtk::ANALYZE_AXIAL:
	orientString = "LAS";
	break;
      case cmtk::ANALYZE_AXIAL_FLIP:
	orientString = "LPS";
	break;
      case cmtk::ANALYZE_CORONAL:
	orientString = "LSA";
	break;
      case cmtk::ANALYZE_CORONAL_FLIP:
	orientString = "LIA";
	break;
      case cmtk::ANALYZE_SAGITTAL:
	orientString = "ASL";
	break;
      case cmtk::ANALYZE_SAGITTAL_FLIP:
	orientString = "AIL";
	break;
      }
    }
  
  UniformVolume::SmartPtr volume( new UniformVolume( dims, UniformVolume::CoordinateVectorType( size ) ) );
  volume->SetMetaInfo( META_IMAGE_ORIENTATION, orientString );
  volume->SetMetaInfo( META_IMAGE_ORIENTATION_ORIGINAL, orientString );

  // Analyze is medical data, which we always treat in RAS space.
  volume->SetMetaInfo( META_SPACE, orientString );
  volume->SetMetaInfo( META_SPACE_ORIGINAL, orientString );
  volume->ChangeCoordinateSpace( AnatomicalOrientation::ORIENTATION_STANDARD );

  // fill "description" header field
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
    case cmtk::ANALYZE_TYPE_NONE:
    case cmtk::ANALYZE_TYPE_BINARY:
    case cmtk::ANALYZE_TYPE_COMPLEX:
    case cmtk::ANALYZE_TYPE_RGB:
    case cmtk::ANALYZE_TYPE_ALL:
    default:
      StdErr.printf( "ERROR: unsupported data type %d in Analyze file %s\n", header.GetField<short>( 70 ), pathHdr.c_str() );
      return volume;
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
    case cmtk::ANALYZE_TYPE_USHORT:
      dtype = TYPE_USHORT;
      break;
    case cmtk::ANALYZE_TYPE_UINT:
      dtype = TYPE_UINT;
      break;
    }  
  
  size_t offset = static_cast<size_t>( header.GetField<float>( 108 ) );
  
  std::string pathImg = pathHdr;
  const size_t period = pathImg.rfind( ".hdr" );
  if ( period != std::string::npos ) 
    pathImg.replace( period, 4, ".img" );
  
  CompressedStream stream( pathImg );
  if ( stream.IsValid() ) 
    {
    stream.Seek( offset, SEEK_CUR );
    
    TypedArray::SmartPtr data( TypedArray::Create( dtype, volume->GetNumberOfPixels() ) );
    if ( data->GetDataSize() == stream.Read( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize() ) )
      {
#ifdef WORDS_BIGENDIAN
      if ( ! bigEndian ) data->ChangeEndianness();
#else
      if ( bigEndian ) data->ChangeEndianness();
#endif
      
      volume->SetData( data );
      }
    else
      {
      StdErr << "ERROR: could not read " << data->GetDataSize() << " pixels from Analyze image file " << pathImg << "\n";
      }
    } 
  else
    {
    StdErr << "ERROR: could not open Analyze image file " << pathImg << "\n";
    }
  
  return volume;
}

void
VolumeFromFile::WriteAnalyzeHdr
( const std::string& pathHdr, const UniformVolume& volume )
{
  UniformVolume::SmartPtr writeVolume( volume.Clone() );
  if ( writeVolume->MetaKeyExists( META_SPACE_ORIGINAL ) )
    writeVolume->ChangeCoordinateSpace( writeVolume->GetMetaInfo( META_SPACE_ORIGINAL ) );

  std::string currentOrientation = writeVolume->GetMetaInfo( META_IMAGE_ORIENTATION );
  if ( currentOrientation == "" )
    {
    currentOrientation = "LAS"; // default: write as is, axial tag, no reorientation.
    }
  std::string originalOrientation = writeVolume->GetMetaInfo( META_IMAGE_ORIENTATION_ORIGINAL );
  if ( originalOrientation == "" )
    {
    originalOrientation = currentOrientation;
    }
  
  // try to write something as close as possible to original orientation
  const char *const supportedOrientations[] = { "LAS", "LSA", "ASL", NULL };
  const char* writeOrientation = AnatomicalOrientation::GetClosestOrientation( originalOrientation.c_str(), supportedOrientations );

  if ( getenv( CMTK_LEGACY_ANALYZE_IO ) || getenv( IGS_LEGACY_ANALYZE_IO ) )
    {
    const char *const supportedOrientationsLegacy[] = { "RAS", "RIP", "AIR", NULL };
    writeOrientation = AnatomicalOrientation::GetClosestOrientation( originalOrientation.c_str(), supportedOrientationsLegacy );
    }
  
  UniformVolume::SmartPtr reorientedVolume;
  if ( strcmp( writeOrientation, currentOrientation.c_str() ) )
    {
    DebugOutput( 2 ) << "INFO: WriteAnalyzeHdr will reorient output volume from '" << currentOrientation << "' to '" << writeOrientation << "'\n";
    reorientedVolume = UniformVolume::SmartPtr( volume.GetReoriented( writeOrientation ) );
    writeVolume = reorientedVolume;
    }
  
  const TypedArray* data = writeVolume->GetData();
  if ( ! data ) return;

  char buffer[348];
  memset( buffer, 0, sizeof( buffer ) );
#ifdef WORDS_BIGENDIAN
  FileHeader header( buffer, true /*bigEndian*/ );
#else
  FileHeader header( buffer, false /*bigEndian*/ );
#endif

  header.StoreField<int>(0, 348 ); // header size
  header.StoreField<int>( 32, 16384 ); // extents
  header.StoreField<short>( 36, 0 ); // session error
  header.StoreField<char>( 38, 'r' ); // regular

  // ndims
  header.StoreField<short>( 40, 4 );

  // dimensions
  header.StoreField<short>( 42, writeVolume->GetDims()[AXIS_X] );
  header.StoreField<short>( 44, writeVolume->GetDims()[AXIS_Y] );
  header.StoreField<short>( 46, writeVolume->GetDims()[AXIS_Z] );
  header.StoreField<short>( 48, 1 ); // write dims 3-7
  header.StoreField<short>( 50, 0 ); // just for safety
  header.StoreField<short>( 52, 0 ); // just for safety
  header.StoreField<short>( 54, 0 ); // just for safety
  header.StoreField<short>( 56, 0 ); // just for safety

  header.StoreField<float>( 68, 0.0 ); // vox_offset
  switch ( data->GetType() ) 
    {
    default:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_NONE );
      header.StoreField<short>( 72, 0 );
    case TYPE_BYTE:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_UNSIGNED_CHAR );
      header.StoreField<short>( 72, 8 * sizeof( byte ) );
      break;
    case TYPE_SHORT:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_SIGNED_SHORT );
      header.StoreField<short>( 72, 8 * sizeof( short ) );
      break;
    case TYPE_USHORT:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_USHORT );
      header.StoreField<short>( 72, 8 * sizeof( unsigned short ) );
      break;
    case TYPE_INT:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_SIGNED_INT );
      header.StoreField<short>( 72, 8 * sizeof( signed int ) );
      break;
    case TYPE_UINT:
      header.StoreField<short>( 70, cmtk::ANALYZE_TYPE_UINT );
      header.StoreField<short>( 72, 8 * sizeof( unsigned int ) );
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
  
  header.StoreField<float>( 80, (float)writeVolume->m_Delta[AXIS_X] );
  header.StoreField<float>( 84, (float)writeVolume->m_Delta[AXIS_Y] );
  header.StoreField<float>( 88, (float)writeVolume->m_Delta[AXIS_Z] );
  header.StoreField<float>( 92, 1.0f ); // write sizes in dims 3 and
  header.StoreField<float>( 96, 1.0f ); // 4 just to be safe

  // set zero offset for binary file.
  header.StoreField<float>( 108, 0.0f );

  // determine data range;
  const Types::DataItemRange dataRange = data->GetRange();

  header.StoreField<float>( 124, static_cast<float>( dataRange.m_UpperBound ) ); // cal_max
  header.StoreField<float>( 128, static_cast<float>( dataRange.m_LowerBound ) ); // cal_min

  header.StoreField<int>( 140, static_cast<int>( dataRange.m_UpperBound ) );
  header.StoreField<int>( 144, static_cast<int>( dataRange.m_LowerBound ) );

  if ( volume.MetaKeyExists( META_IMAGE_DESCRIPTION ) )
    header.StoreFieldString( 148, volume.GetMetaInfo( META_IMAGE_DESCRIPTION ).c_str(), 80 );
  
  if ( getenv( CMTK_LEGACY_ANALYZE_IO ) || getenv( IGS_LEGACY_ANALYZE_IO ) )
    {
    // slice orientation always axial from caudal to cranial
    header.StoreField<byte>( 252, cmtk::ANALYZE_AXIAL );
    header.StoreField<byte>( 254, 0 ); //set Nifti sform code to 0.
    }
  else
    {
    if ( !strcmp( writeOrientation, "LAS" ) )
      header.StoreField<byte>( 252, cmtk::ANALYZE_AXIAL );
    else if ( !strcmp( writeOrientation, "LSA" ) )
      header.StoreField<byte>( 252, cmtk::ANALYZE_CORONAL );
    else if ( !strcmp( writeOrientation, "ASL" ) )
      header.StoreField<byte>( 252, cmtk::ANALYZE_SAGITTAL );
    header.StoreField<byte>( 254, 0 ); //set Nifti sform code to 0.

    // mark this as "new" SRI Analyze image.
    header.StoreFieldString( 344, "SRI1", 4 );
    }

  // write binary data
  std::string pathImg = pathHdr;
  const size_t period = pathImg.rfind( ".hdr" );
  if ( period != std::string::npos ) 
    pathImg.replace( period, 4, ".img" );
  
  if ( VolumeIO::GetWriteCompressed() )
    {
    struct stat buf;
    if ( ! stat( pathImg.c_str(), &buf ) )
      {
      StdErr << "WARNING: Analyze img file '" << pathImg << "' will be written compressed, but uncompressed file exists!\n";
      }
    
#ifdef _MSC_VER
    const char *const modestr = "w9b";
#else
    const char *const modestr = "w9";
#endif
    
    gzFile imgFile = gzopen( (pathImg + ".gz").c_str(), modestr );
    if ( imgFile ) 
      {
      const size_t dataSize = data->GetItemSize() * data->GetDataSize();
      if ( dataSize != CompressedStream::Zlib::StaticSafeWrite( imgFile, data->GetDataPtr(), dataSize ) )
	{
	StdErr << "WARNING: gzwrite() returned error when writing to " << pathImg << "\n";
	}
      gzclose( imgFile );
      }
    }
  else
    {
#ifdef _MSC_VER
    const char *const modestr = "wb";
#else
    const char *const modestr = "w";
#endif

    FILE *imgFile = fopen( pathImg.c_str(), modestr );
    if ( imgFile ) 
      {
      fwrite( data->GetDataPtr(), data->GetItemSize(), data->GetDataSize(), imgFile );
      fclose( imgFile );
      }
    else
      {
      StdErr << "ERROR: could not open file '" << pathImg << "' for writing\n";
      }
    }

  // write header info
#ifdef _MSC_VER
  FILE *hdrFile = fopen( pathHdr.c_str(), "wb" );
#else
  FILE *hdrFile = fopen( pathHdr.c_str(), "w" );
#endif
  if ( hdrFile ) 
    {
    if ( 348 != fwrite( buffer, 1, 348, hdrFile ) ) 
      {
      StdErr << "ERROR: could not write 348 bytes to header file " << pathHdr << "\n";
      }
    fclose( hdrFile );
    }
  else
    {
    StdErr << "ERROR: could not open file '" << pathHdr << "' for writing\n";
    }
}

} // namespace cmtk
