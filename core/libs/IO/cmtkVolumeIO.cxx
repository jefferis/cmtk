/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
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

#include "cmtkVolumeIO.h"

#include <System/cmtkStrUtility.h>
#include <System/cmtkFileUtils.h>
#include <System/cmtkProgress.h>
#include <System/cmtkMountPoints.h>

#include <IO/cmtkStudy.h>
#include <IO/cmtkClassStream.h>
#include <IO/cmtkVolumeFromFile.h>
#include <IO/cmtkPGM.h>

#include <Base/cmtkTypes.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef HAVE_LIBGEN_H
#  include <libgen.h>
#endif

#include <limits.h>
#include <math.h>
#include <stdio.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Environment variable that turns off writing output images in the array order of the input images.
const char* const CMTK_LEGACY_WRITE_IMAGES_RAS = "CMTK_LEGACY_WRITE_IMAGES_RAS";

/// Static global flag: when true, files are written compressed whenever possible.
bool VolumeIO::WriteCompressedOn = true;

UniformVolume::SmartPtr
VolumeIO::Read( const char* path, const bool verbose )
{
  UniformVolume::SmartPtr volume( NULL );

  const char *translatedPath = MountPoints::Translate( path );

  FileFormatID formatID = FileFormat::Identify( translatedPath );
  switch ( formatID ) 
    {
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
    case FILEFORMAT_NRRD:
      volume = VolumeFromFile::ReadNRRD( translatedPath );
      break;
    default: 
    {
    }
    }
  
  if ( volume ) 
    {
    volume->m_MetaInformation[META_FILEFORMAT_ORIGINAL] = FileFormat::Describe( formatID );
    // for float and double data, automatically recognize Inf as Null Data.
    TypedArray::SmartPtr dataArray = volume->GetData();
    if ( dataArray ) 
      {
      if ( dataArray->GetType() == TYPE_FLOAT )
	{
	const float fInf = MathUtil::GetFloatInf();
	dataArray->SetPaddingPtr( &fInf );
	}
      else
	{
	if ( dataArray->GetType() == TYPE_DOUBLE )
	  {
	  const double dInf = MathUtil::GetDoubleInf();
	  dataArray->SetPaddingPtr( &dInf );
	  }
	}

      for ( size_t i = 0; i < dataArray->GetDataSize(); ++i )
	{
	Types::DataItem v;
	if ( dataArray->Get( v, i ) && MathUtil::IsNaN( v ) )
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
      const Types::DataItemRange range = dataArray->GetRange();
      StdErr.printf( "Data type %s, range [%f .. %f]\n", DataTypeName[ dataArray->GetType() ],
		     static_cast<float>( range.m_LowerBound ), static_cast<float>( range.m_UpperBound ) );
      } 
    else
      {
      StdErr << "Image does not contain valid data.\n";
      }
    }

  if ( volume )
    {
    volume->m_MetaInformation[META_FS_PATH] = path;
    }
  
  return volume;
}

UniformVolume::SmartPtr
VolumeIO::ReadGrid( const char* path, const bool verbose )
{
  UniformVolume::SmartPtr volume( NULL );

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
    volume->m_MetaInformation[META_FS_PATH] = path;
    }
  
  return volume;
}

UniformVolume::SmartPtr
VolumeIO
::ReadGridOriented( const char *path, const char* orientation, const bool verbose )
{
  UniformVolume::SmartPtr volume( Self::ReadGrid( path, verbose ) );
  if ( !volume ) 
    return volume;
  
  const std::string orientationOriginal = volume->m_MetaInformation[META_IMAGE_ORIENTATION];
  if ( orientationOriginal == "" )
    {
    StdErr << "WARNING: image does not have valid orientation meta information; cannot reorient.\n";
    return volume;
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

  return volume;
}

UniformVolume::SmartPtr
VolumeIO
::ReadOriented( const char *path, const char* orientation, const bool verbose )
{
  UniformVolume::SmartPtr volume( VolumeIO::Read( path, verbose ) );
  if ( !volume ) 
    return volume;

  const std::string orientationOriginal = volume->m_MetaInformation[META_IMAGE_ORIENTATION];
  if ( orientationOriginal == "" )
    {
    StdErr << "WARNING: image does not have valid orientation meta information; cannot reorient.\n";
    return volume;
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
  return volume;
}

void 
VolumeIO::Write
( const UniformVolume& volume, const char *pathAndFormat, const bool verbose )
{
  const char* actualPath = pathAndFormat;
  FileFormatID fileFormat = FILEFORMAT_UNKNOWN;

  const char* suffix = strrchr( pathAndFormat, '.' );  
  if ( suffix )
    {
    // check whether we have a compression-related suffix
    if ( ! strcmp( suffix, ".gz" ) )
      {
      // include actual suffix
      while ( suffix != pathAndFormat )
	{
	--suffix;
	if ( *suffix == '.' )
	  break;
	}
      }

    if ( ! strcmp( ".hdr", suffix ) )
      {
      fileFormat = FILEFORMAT_ANALYZE_HDR;
      }
    else
      {
      if ( ! strcmp( ".img", suffix ) || ! strcmp( ".img.gz", suffix ) )
	{
	fileFormat = FILEFORMAT_NIFTI_DETACHED;
	}
      else
	{
	if ( ! strcmp( ".nii", suffix ) || ! strcmp( ".nii.gz", suffix ) )
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
    
    if ( ! strncmp( "ANALYZE", pathAndFormat, formatLength ) ) 
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
    }
#endif

  if ( fileFormat == FILEFORMAT_UNKNOWN )
    {
    StdErr << "Fileformat not recognized; writing single-file NIFTI instead.\n";
    fileFormat = FILEFORMAT_NIFTI_SINGLEFILE;
    }
  
  char absolutePath[PATH_MAX];
  FileUtils::GetAbsolutePath( absolutePath, actualPath );
  
  Write( volume, fileFormat, absolutePath, verbose );
}

void
VolumeIO::Write
( const UniformVolume& volume, const FileFormatID format, const char* path, const bool verbose )
{
  const TypedArray *data = volume.GetData();
  if ( data == NULL ) return;

  const UniformVolume* actualVolume = &volume;

  // if volume was reoriented from its original array order, temporarily reorient back and set actual volume to temporary volume.
  cmtk::UniformVolume::SmartConstPtr reorientedVolume;

  if ( !getenv( CMTK_LEGACY_WRITE_IMAGES_RAS ) )
    {
    if ( volume.MetaKeyExists( cmtk::META_IMAGE_ORIENTATION_ORIGINAL ) &&
	 (volume.GetMetaInfo( cmtk::META_IMAGE_ORIENTATION ) != volume.GetMetaInfo( cmtk::META_IMAGE_ORIENTATION_ORIGINAL ) ) )
      {
      reorientedVolume = cmtk::UniformVolume::SmartConstPtr( volume.GetReoriented( volume.GetMetaInfo( cmtk::META_IMAGE_ORIENTATION_ORIGINAL ).c_str() ) );
      actualVolume = reorientedVolume;
      }
    }
    
  switch ( format ) 
    {
    case FILEFORMAT_ANALYZE_HDR: 
    {
    VolumeFromFile::WriteAnalyzeHdr( path, *actualVolume, verbose );
    break;
    }
    case FILEFORMAT_NIFTI_DETACHED: 
    case FILEFORMAT_NIFTI_SINGLEFILE: 
    {
    VolumeFromFile::WriteNifti( path, *actualVolume, verbose );
    break;
    }
    case FILEFORMAT_METAIMAGE: 
    {
    VolumeFromFile::WriteMetaImage( path, *actualVolume );
    break;
    }
    case FILEFORMAT_NRRD: 
    {
    VolumeFromFile::WriteNRRD( path, *actualVolume, verbose );
    break;
    }
    break;
    default:
      break;
    }
  
  volume.m_MetaInformation[META_FS_PATH] = path;
}

VolumeIO::Initializer::Initializer()
{
  if ( getenv( "IGS_WRITE_UNCOMPRESSED" ) || getenv( "CMTK_WRITE_UNCOMPRESSED" ) )
    VolumeIO::SetWriteCompressedOff();
}

VolumeIO::Initializer VolumeIO::Initializer::Instance;

} // namespace cmtk
