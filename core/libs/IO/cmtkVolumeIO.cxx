/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#include <cmtkVolumeIO.h>

#include <cmtkMountPoints.h>
#include <cmtkStudy.h>
#include <cmtkClassStream.h>
#include <cmtkVolumeFromStudy.h>
#include <cmtkVolumeFromFile.h>
#include <cmtkPGM.h>

#include <cmtkStrUtility.h>
#include <cmtkFileUtil.h>
#include <cmtkProgress.h>
#include <cmtkTypes.h>

#include <stdio.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef HAVE_LIBGEN_H
#  include <libgen.h>
#endif

#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif

#include <math.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

bool VolumeIO::WriteCompressedOn = true;

UniformVolume* 
VolumeIO::Read( const char* path, const bool verbose )
{
  UniformVolume *volume = NULL;

  const char *translatedPath = MountPoints::Translate( path );

  FileFormatID formatID = FileFormat::Identify( translatedPath );
  switch ( formatID ) 
    {
    case FILEFORMAT_STUDY: 
    {
    Study *study = NULL;
    ClassStream stream( translatedPath, "images", ClassStream::READ );
    stream >> study;
    if ( study ) 
      {
      volume = VolumeFromStudy::Read( study, verbose );
      delete study;
      }
    }
    break;
    case FILEFORMAT_DICOM: // (hopefully) multi-slice DICOM
      volume = VolumeFromFile::ReadDICOM( translatedPath );
      break;
    case FILEFORMAT_VANDERBILT:
      volume = VolumeFromFile::ReadVanderbilt( translatedPath );
      break;
    case FILEFORMAT_BIORAD:
      volume = VolumeFromFile::ReadBioRad( translatedPath );
      break;
    case FILEFORMAT_ANALYZE_HDR:
      volume = VolumeFromFile::ReadAnalyzeHdr( translatedPath, false /*bigendian*/, true /*readData*/ );
      break;
    case FILEFORMAT_ANALYZE_HDR_BIGENDIAN:
      volume = VolumeFromFile::ReadAnalyzeHdr( translatedPath, true /*bigendian*/, true /*readData*/ );
      break;
    case FILEFORMAT_NIFTI_SINGLEFILE:
      volume = VolumeFromFile::ReadNifti( translatedPath, false /*detached*/, true /*readData*/ );
      break;
    case FILEFORMAT_NIFTI_DETACHED:
      volume = VolumeFromFile::ReadNifti( translatedPath, true /*detached*/, true /*readData*/ );
      break;
    case FILEFORMAT_ANALYZE_AVW:
      volume = VolumeFromFile::ReadAnalyzeAVW( translatedPath );
      break;
    case FILEFORMAT_NRRD:
      volume = VolumeFromFile::ReadNRRD( translatedPath );
      break;
    default: 
    {
    }
    }
  
  if ( volume ) 
    {
    volume->m_MetaInformation[CMTK_META_FILEFORMAT_ORIGINAL] = FileFormat::Describe( formatID );
    // for float and double data, automatically recognize Inf as Null Data.
    TypedArray::SmartPtr dataArray = volume->GetData();
    if ( !dataArray.IsNull() ) 
      {
      if ( dataArray->GetType() == TYPE_FLOAT )
	dataArray->SetPaddingPtr( &CMTK_FLOAT_INF );
      else
	if ( dataArray->GetType() == TYPE_DOUBLE )
	  dataArray->SetPaddingPtr( &CMTK_DOUBLE_INF );

      for ( size_t i = 0; i < dataArray->GetDataSize(); ++i )
	{
	Types::DataItem v;
	if ( dataArray->Get( v, i ) && isnan( v ) )
	  dataArray->SetPaddingAt( i );
	}
      }
    }
  
  if ( verbose && volume ) 
    {
    StdErr.printf( "%s\nRead %d x %d x %d voxels [%f x %f x %f mm total size].\n", path,
		      volume->GetDims()[0], volume->GetDims()[1], volume->GetDims()[2],
		      volume->Size[0], volume->Size[1], volume->Size[2] );
    
    const TypedArray* dataArray = volume->GetData();
    if ( dataArray ) 
      {
      Types::DataItem min = 0, max = 0;
      if ( dataArray->GetRange( min, max ) ) 
	{
	StdErr.printf( "Data type %s, range [%f .. %f]\n", DataTypeName[ dataArray->GetType() ],
			  static_cast<float>( min ), static_cast<float>( max ) );
	} 
      else
	{
	StdErr << "Image does not contain valid data.\n";
	}
      }
    }

  if ( volume )
    {
    volume->m_MetaInformation[CMTK_META_FS_PATH] = path;
    }
  
  return volume;
}

UniformVolume* 
VolumeIO::ReadGrid( const char* path, const bool verbose )
{
  UniformVolume *volume = NULL;

  const char *translatedPath = MountPoints::Translate( path );

  switch ( FileFormat::Identify( translatedPath ) ) 
    {
    case FILEFORMAT_ANALYZE_HDR:
      volume = VolumeFromFile::ReadAnalyzeHdr( translatedPath, false /*bigendian*/, false /*readData*/ );
      break;
    case FILEFORMAT_ANALYZE_HDR_BIGENDIAN:
      volume = VolumeFromFile::ReadAnalyzeHdr( translatedPath, true /*bigendian*/, false /*readData*/ );
      break;
    case FILEFORMAT_NIFTI_SINGLEFILE:
      volume = VolumeFromFile::ReadNifti( translatedPath, false /*detached*/, false /*readData*/ );
      break;
    case FILEFORMAT_NIFTI_DETACHED:
      volume = VolumeFromFile::ReadNifti( translatedPath, true /*detached*/, false /*readData*/ );
      break;
    default: 
    {
    // For now, default to full reader for other image file formats
    volume = VolumeIO::Read( path, verbose );
    }
    }
  
  if ( verbose && volume ) 
    {
    StdErr.printf( "%s\nRead %d x %d x %d voxels [%f x %f x %f mm total size].\n", path,
		      volume->GetDims()[0], volume->GetDims()[1], volume->GetDims()[2],
		      volume->Size[0], volume->Size[1], volume->Size[2] );
    }
  
  if ( volume )
    {
    volume->m_MetaInformation[CMTK_META_FS_PATH] = path;
    }
  
  return volume;
}

UniformVolume* 
VolumeIO
::ReadGridOriented( const char *path, const char* orientation, const bool verbose )
{
  UniformVolume::SmartPtr volume( Self::ReadGrid( path, verbose ) );
  if ( !volume ) return NULL;
  
  const std::string orientationOriginal = volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
  if ( orientationOriginal == "" )
    {
    StdErr << "WARNING: image does not have valid orientation meta information; cannot reorient.\n";
    return volume.ReleasePtr();
    }
  else
    {
    if ( orientationOriginal != orientation )
      {
      if ( verbose )
	{
	StdErr << "INFO: reorienting image from '" << orientationOriginal << "' to '" << orientation << "'\n";
	}
      
      return volume->GetReoriented( orientation );
      }
    }

  return volume.ReleasePtr();
}

UniformVolume* 
VolumeIO
::ReadOriented( const char *path, const char* orientation, const bool verbose )
{
  UniformVolume::SmartPtr volume( VolumeIO::Read( path, verbose ) );
  if ( !volume ) return NULL;

  const std::string orientationOriginal = volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
  if ( orientationOriginal == "" )
    {
    StdErr << "WARNING: image does not have valid orientation meta information; cannot reorient.\n";
    return volume.ReleasePtr();
    }
  else
    {
    if ( orientationOriginal != orientation )
      {
      if ( verbose )
	{
	StdErr << "INFO: reorienting image from '" << orientationOriginal << "' to '" << orientation << "'\n";
	}
      
      return volume->GetReoriented( orientation );
      }
    }
  return volume.ReleasePtr();
}

void 
VolumeIO::Write
( const UniformVolume* volume, const char *pathAndFormat, const bool verbose )
{
  const char* actualPath = pathAndFormat;
  FileFormatID fileFormat = FILEFORMAT_UNKNOWN;

  const char* suffix = strrchr( pathAndFormat, '.' );
  if ( suffix )
    {
    if ( ! strcmp( ".hdr", suffix ) )
      {
      fileFormat = FILEFORMAT_ANALYZE_HDR;
      }
    else
      {
      if ( ! strcmp( ".img", suffix ) )
	{
	fileFormat = FILEFORMAT_NIFTI_DETACHED;
	}
      else
	{
	if ( ! strcmp( ".nii", suffix ) )
	  {
	  fileFormat = FILEFORMAT_NIFTI_SINGLEFILE;
	  }
	else
	  {
	  if ( ! strcmp( ".mha", suffix ) )
	    {
	    fileFormat = FILEFORMAT_METAIMAGE;
	    }
	  else 
	    {
	    if ( ! strcmp( ".nrrd", suffix ) || ! strcmp( ".nhdr", suffix ) )
	      {
	      fileFormat = FILEFORMAT_NRRD;
	      }
	    }
	  }
	}
      }
    }
 
#ifndef _MSC_VER
  const char* colon = strchr( pathAndFormat, ':' );
  if ( colon != NULL ) 
    {
    actualPath = colon+1;
    unsigned int formatLength = colon - pathAndFormat - 1;
    
    if ( ! strncmp( "RAW3D", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_RAW3D;
      } 
    else if ( ! strncmp( "ANALYZE", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_ANALYZE_HDR;
      } 
    else if ( ! strncmp( "NIFTI", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_NIFTI_SINGLEFILE;
      } 
    else if ( ! strncmp( "NRRD", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_NRRD;
      } 
    else if ( ! strncmp( "METAIMAGE", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_METAIMAGE;
      } 
    else if ( ! strncmp( "NRRD", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_NRRD;
      } 
    else if ( ! strncmp( "VANDERBILT", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_VANDERBILT;
      } 
    else if ( ! strncmp( "RAW", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_RAW;
      } 
    else if ( ! strncmp( "PGM", pathAndFormat, formatLength ) ) 
      {
      fileFormat = FILEFORMAT_PGM;
      }
    }
#endif

  if ( fileFormat == FILEFORMAT_UNKNOWN )
    {
    StdErr << "Fileformat not recognized; writing RAW3D instead.\n";
    fileFormat = FILEFORMAT_RAW3D;
    }
  
  char absolutePath[PATH_MAX];
  FileUtils::GetAbsolutePath( absolutePath, actualPath );
  
  char studyPath[PATH_MAX];
  strcat( strcpy( studyPath, StrDir( absolutePath ) ), ".study" );

  Write( volume, fileFormat, absolutePath, studyPath, verbose );
}

void
VolumeIO::Write
( const UniformVolume* volume, const FileFormatID format, const char* path, const char* studyPath, const bool verbose )
{
  if ( volume == NULL ) return;
 
  const TypedArray *data = volume->GetData();
  if ( data == NULL ) return;

  int planeSize = volume->m_Dims[0] * volume->m_Dims[1];
  ScalarImage image( volume->m_Dims[0], volume->m_Dims[1] );
  image.SetPixelSize( volume->m_Delta[AXIS_X], volume->m_Delta[AXIS_Y] );

  FileUtils::RecursiveMkPrefixDir( path );
  char *dirName = strdup( StrDir( path ) );
  char *baseName = strdup( StrFName( path ) );

  switch ( format ) 
    {
    case FILEFORMAT_PGM: 
    {
    int offset = 0;
    TypedStream *studyStream = NULL;
    if ( studyPath ) 
      {
      studyStream = CreateStudy( studyPath, volume, dirName, FILEFORMAT_PGM, verbose );
      }
    for ( int k = 0; k < volume->GetDims()[2]; ++k, offset += planeSize ) 
      {
      if ( verbose ) 
	{
	fprintf( stderr, "Writing slice %d...\r", k );
	}
      image.SetPixelData( TypedArray::SmartPtr( data->CloneSubArray( offset, planeSize ) ) );
      char* fname = MakeSliceFileName( path, "pgm", k );
      PGM::Write( fname, &image );
      if ( studyStream ) 
	{
	WriteStudyImageEntry( studyStream, StrFName( fname ), "PRIMARY", volume->m_Delta[AXIS_X], volume->m_Delta[AXIS_Y], volume->GetPlaneCoord( AXIS_Z, k ) );
	}
      }
    if ( studyStream ) 
      {
      FinishStudy( studyStream );
      }
    }
    break;
    case FILEFORMAT_VANDERBILT: 
    {
    // get image data as short array
    short *shortData = static_cast<short*>( data->ConvertArray( TYPE_SHORT ) );
#ifndef WORDS_BIGENDIAN
    // change endianness from Sun to whatever we're currently on.
    //    data->ChangeEndianness();
#endif
    WriteData( path, shortData, data->GetDataSize(), sizeof( *shortData ) );
    
    // free using allocator used by the array object
    data->Free( shortData );
    }
    break;
    case FILEFORMAT_ANALYZE_HDR: 
    {
    VolumeFromFile::WriteAnalyzeHdr( path, volume, verbose );
    break;
    }
    case FILEFORMAT_NIFTI_DETACHED: 
    case FILEFORMAT_NIFTI_SINGLEFILE: 
    {
    VolumeFromFile::WriteNifti( path, volume, verbose );
    break;
    }
    case FILEFORMAT_METAIMAGE: 
    {
    VolumeFromFile::WriteMetaImage( path, volume );
    break;
    }
    case FILEFORMAT_NRRD: 
    {
    VolumeFromFile::WriteNRRD( path, volume, verbose );
    break;
    }
    break;
    case FILEFORMAT_RAW3D: 
    {
#ifndef WORDS_BIGENDIAN
    // change endianness from Sun to whatever we're currently on.
    //    data->ChangeEndianness();
#endif
    WriteData( path, data->GetDataPtr(), data->GetDataSize(), 
	       data->GetItemSize() );
    if ( studyPath ) 
      {
      TypedStream *studyStream = CreateStudy( studyPath, volume, dirName, FILEFORMAT_RAW3D, verbose );
      WriteStudyImageEntry( studyStream, baseName, "PRIMARY" );
      FinishStudy( studyStream );
      }
    }
    break;
    default:
      break;
    }
  
  free( dirName );
  free( baseName );
  
  volume->m_MetaInformation[CMTK_META_FS_PATH] = path;
}

bool 
VolumeIO::WriteData
( const char* path, const void *dataPtr, const size_t numberOfItems, const size_t itemSize )
{
  FILE *fp = fopen( path, "w" );
  if ( fp ) 
    {
    Progress::Begin( 0, 1 + numberOfItems / (1<<20), 1, "Writing volume data" );
    unsigned int step = 0;
    
    // if bigger than 1 MB, write in 1 MB chunks to allow progress reports
    const char* charDataPtr = static_cast<const char*>( dataPtr );
    size_t itemsLeft = numberOfItems;
    while ( itemsLeft ) 
      {
      if ( itemsLeft < (1<<20) ) 
	{
	fwrite( charDataPtr, itemSize, itemsLeft, fp );
	itemsLeft = 0;
	} 
      else
	{
	fwrite( charDataPtr, itemSize, (1<<20), fp );
	charDataPtr += (1<<20) * itemSize;
	itemsLeft -= (1<<20);
	}
      Progress::SetProgress( ++step );
      }
    fclose( fp );
    } 
  else
    {
    return false;
    }
  return true;
}

char* 
VolumeIO::MakeSliceFileName
( const char *path, const char* suffix, const int index )
{
  static char fullname[PATH_MAX];
  snprintf( fullname, sizeof( fullname ), path, index, suffix );
  return fullname;
}

TypedStream*
VolumeIO::CreateStudy
( const char* studyPath, const UniformVolume* volume, const char* path, 
  const FileFormatID fileFormat, const bool verbose )
{
  const TypedArray *volumeData = volume->GetData();
  TypedStream* studyStream = new TypedStream( studyPath, "images", TYPEDSTREAM_WRITE );
  
  if ( studyStream && !studyStream->IsValid() ) 
    {
    fprintf( stderr, "Failed to create study %s.\n", studyPath );
    delete studyStream;
    studyStream = NULL;
    }
  
  if ( studyStream )
    {
    if ( verbose ) 
      {
      fprintf( stderr, "Creating study %s.\n", studyPath );
      }
    studyStream->Begin( "imageserie" );
    studyStream->WriteString( "name", "PRIMARY" );
    studyStream->WriteString( "imagepath", path );
    studyStream->WriteString( "iconpath", path );

    switch ( fileFormat ) 
      {
      case FILEFORMAT_DICOM:
	studyStream->WriteString( "format", "DICOM" );
	break;
      case FILEFORMAT_PGM:
	studyStream->WriteString( "format", "PGM" );
	break;
      case FILEFORMAT_RAW:
	studyStream->WriteString( "format", "RAW-DATA" );
	break;
      case FILEFORMAT_RAW3D:
	studyStream->WriteString( "format", "RAW3D" );
	break;
      default:
	break;
      }
    
    studyStream->WriteInt( "width", volume->GetDims( AXIS_X ) );
    studyStream->WriteInt( "height", volume->GetDims( AXIS_Y ) );
    if ( fileFormat == FILEFORMAT_RAW3D )
      studyStream->WriteInt( "depth", volume->GetDims( AXIS_Z ) );
    
    studyStream->WriteString( "direction", "CaudalCranial" );
    studyStream->WriteBool( "custom", 1 );
    studyStream->WriteDouble( "calibrationx", volume->m_Delta[AXIS_X] );
    studyStream->WriteDouble( "calibrationy", volume->m_Delta[AXIS_Y] );
    studyStream->WriteDouble( "slicedistance", volume->m_Delta[AXIS_Z] );

    if ( volumeData ) 
      {
      Types::DataItem min, max;
      volumeData->GetRange( min, max );
      studyStream->WriteInt( "minimum", static_cast<int>( min ) );
      studyStream->WriteInt( "maximum", static_cast<int>( max ) );

      studyStream->WriteInt( "offset", 0 );
      studyStream->WriteInt( "bytesperpixel", volumeData->GetItemSize() );
      studyStream->WriteBool( "signed", (min < 0) );

#ifdef WORDS_BIGENDIAN
      studyStream->WriteBool( "bigendian", true );
      studyStream->WriteBool( "littleendian", false );
      studyStream->WriteBool( "swapbytes", false );
#else
      studyStream->WriteBool( "bigendian", false );
      studyStream->WriteBool( "littleendian", true );
      studyStream->WriteBool( "swapbytes", true );
#endif
      
      switch ( volumeData->GetDataClass() ) 
	{
	case DATACLASS_GREY:
	  studyStream->WriteString( "dataclass", "Grey" );
	  break;
	case DATACLASS_LABEL:
	  studyStream->WriteString( "dataclass", "Label" );
	  break;
	default:
	  studyStream->WriteString( "dataclass", "Unknown" );
	  break;
	}
      
      if ( volumeData->GetPaddingFlag() ) 
	{
	studyStream->WriteBool( "padding", true );
	studyStream->WriteDouble( "padding_value", volumeData->GetPaddingValue() );
	}
      }
    
    studyStream->End( TYPEDSTREAM_FLUSH );
    }
  return studyStream;
}

void
VolumeIO::WriteStudyImageEntry
( TypedStream *const studyStream, const char* fname, const char* serie, 
  const Types::Coordinate calibrationx, const Types::Coordinate calibrationy,
  const Types::Coordinate tablepos )
{
  if ( studyStream->IsValid() ) 
    {    
    studyStream->Begin( "image" );
    studyStream->WriteString( "serie", serie );
    studyStream->WriteString( "name", fname );
    studyStream->WriteDouble( "calibrationx", calibrationx );
    studyStream->WriteDouble( "calibrationy", calibrationy );
    studyStream->WriteDouble( "tablepos", tablepos );
    studyStream->End( TYPEDSTREAM_FLUSH );
    }
}

void VolumeIO::FinishStudy( TypedStream *const studyStream )
{
  studyStream->Close();
  delete studyStream;
}

VolumeIO::Initializer::Initializer()
{
  if ( getenv( "IGS_WRITE_UNCOMPRESSED" ) || getenv( "CMTK_WRITE_UNCOMPRESSED" ) )
    VolumeIO::SetWriteCompressedOff();
}

VolumeIO::Initializer VolumeIO::Initializer::Instance;

} // namespace cmtk
