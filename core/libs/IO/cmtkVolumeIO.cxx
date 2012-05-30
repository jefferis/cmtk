/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkVolumeIO.h"

#include <System/cmtkExitException.h>
#include <System/cmtkStrUtility.h>
#include <System/cmtkFileUtils.h>
#include <System/cmtkProgress.h>
#include <System/cmtkMountPoints.h>
#include <System/cmtkDebugOutput.h>

#include <IO/cmtkStudy.h>
#include <IO/cmtkClassStreamInput.h>
#include <IO/cmtkClassStreamOutput.h>
#include <IO/cmtkVolumeFromFile.h>

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

#include <string>

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
VolumeIO::Read( const std::string& path )
{
  UniformVolume::SmartPtr volume( NULL );

  const std::string translatedPath = MountPoints::Translate( path );

  const FileFormatID formatID = FileFormat::Identify( translatedPath );
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
    default: {}
    }
  
  if ( ! volume )
    {
    StdErr << "ERROR: could not read image geometry from " << path << "\n";
    throw ExitException( 1 );
    }
  
  volume->SetMetaInfo( META_FS_PATH, path );
  volume->SetMetaInfo( META_FILEFORMAT_ORIGINAL, FileFormat::Describe( formatID ) );
  DebugOutput( 3 ).GetStream().printf( "%s\nRead %d x %d x %d voxels [%f x %f x %f mm total size].\n", path.c_str(), 
				       volume->GetDims()[0], volume->GetDims()[1], volume->GetDims()[2], volume->m_Size[0], volume->m_Size[1], volume->m_Size[2] );
  
  const TypedArray* dataArray = volume->GetData();
  if ( ! dataArray )
    {
    StdErr << "ERROR: could not read image data from " << path << "\n";
    throw ExitException( 1 );
    }
  
  const Types::DataItemRange range = dataArray->GetRange();
  DebugOutput( 3 ).GetStream().printf( "Data type %s, range [%f .. %f]\n", DataTypeName[ dataArray->GetType() ], static_cast<float>( range.m_LowerBound ), static_cast<float>( range.m_UpperBound ) );

  return volume;
}

UniformVolume::SmartPtr
VolumeIO::ReadGrid( const std::string& path )
{
  UniformVolume::SmartPtr volume( NULL );

  const std::string translatedPath = MountPoints::Translate( path );

  const FileFormatID formatID = FileFormat::Identify( translatedPath );
  try
    {
    switch ( formatID ) 
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
      default: {
      // For now, default to full reader for other image file formats
      volume = VolumeIO::Read( path );
      }
      }
    }
  catch (...) {}

  if ( ! volume )
    {
    StdErr << "ERROR: could not read image from " << path << "\n";
    throw ExitException( 1 );
    }
  
  DebugOutput( 3 ).GetStream().printf( "%s\nRead %d x %d x %d voxels [%f x %f x %f mm total size].\n", path.c_str(), 
				       volume->GetDims()[0], volume->GetDims()[1], volume->GetDims()[2], volume->m_Size[0], volume->m_Size[1], volume->m_Size[2] );

  volume->SetMetaInfo( META_FS_PATH, path );
  volume->SetMetaInfo( META_FILEFORMAT_ORIGINAL, FileFormat::Describe( formatID ) );
  
  return volume;
}

UniformVolume::SmartPtr
VolumeIO
::ReadGridOriented( const std::string& path, const char* orientation )
{
  UniformVolume::SmartPtr volume( Self::ReadGrid( path ) );
  
  const std::string orientationOriginal = volume->GetMetaInfo( META_IMAGE_ORIENTATION );
  if ( orientationOriginal == "" )
    {
    StdErr << "WARNING: image does not have valid orientation meta information; cannot reorient.\n";
    return volume;
    }
  else
    {
    if ( orientationOriginal != orientation )
      {
      DebugOutput( 3 ) << "Reorienting image from '" << orientationOriginal << "' to '" << orientation << "'\n";
      return volume->GetReoriented( orientation );
      }
    }

  return volume;
}

UniformVolume::SmartPtr
VolumeIO
::ReadOriented( const char *path, const char* orientation )
{
  UniformVolume::SmartPtr volume( VolumeIO::Read( path ) );

  const std::string orientationOriginal = volume->GetMetaInfo( META_IMAGE_ORIENTATION );
  if ( orientationOriginal == "" )
    {
    StdErr << "WARNING: image does not have valid orientation meta information; cannot reorient.\n";
    return volume;
    }
  else
    {
    if ( orientationOriginal != orientation )
      {
      DebugOutput( 3 ) << "INFO: reorienting image from '" << orientationOriginal << "' to '" << orientation << "'\n";      
      return volume->GetReoriented( orientation );
      }
    }
  return volume;
}

void 
VolumeIO::Write
( const UniformVolume& volume, const std::string& pathAndFormat )
{
  std::string actualPath = pathAndFormat;
  FileFormatID fileFormat = FILEFORMAT_UNKNOWN;

  const size_t period = pathAndFormat.rfind( '.' );
  if ( period != std::string::npos )
    {
    std::string suffix = pathAndFormat.substr( period );  

    // check whether we have a compression-related suffix
    if ( suffix == ".gz" )
      {
      // include actual suffix
      const size_t period2 = pathAndFormat.rfind( '.', period-1 );
      suffix = pathAndFormat.substr( period2, period-period2 );
      }

    if ( suffix == ".hdr" )
      {
      fileFormat = FILEFORMAT_ANALYZE_HDR;
      }
    else
      {
      if ( suffix == ".img" )
	{
	fileFormat = FILEFORMAT_NIFTI_DETACHED;
	}
      else
	{
	if ( suffix == ".nii" )
	  {
	  fileFormat = FILEFORMAT_NIFTI_SINGLEFILE;
	  }
	else
	  {
	  if ( suffix == ".mha" )
	    {
	    fileFormat = FILEFORMAT_METAIMAGE;
	    }
	  else 
	    {
	    if ( ( suffix == ".nrrd") || (suffix == ".nhdr") )
	      {
	      fileFormat = FILEFORMAT_NRRD;
	      }
	    }
	  }
	}
      }
    }
  
#ifndef _MSC_VER
  const size_t colon = pathAndFormat.find( ':' );
  if ( colon != std::string::npos ) 
    {
    actualPath = pathAndFormat.substr( colon+1 );
    const std::string format = pathAndFormat.substr( colon-1 );
    
    if ( format == "ANALYZE" ) 
      {
      fileFormat = FILEFORMAT_ANALYZE_HDR;
      } 
    else if ( format == "NIFTI" ) 
      {
      fileFormat = FILEFORMAT_NIFTI_SINGLEFILE;
      } 
    else if ( format == "NRRD" ) 
      {
      fileFormat = FILEFORMAT_NRRD;
      } 
    else if ( format == "METAIMAGE" ) 
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
  
  const std::string absolutePath = FileUtils::GetAbsolutePath( actualPath );
  Write( volume, fileFormat, absolutePath );
}

void
VolumeIO::Write
( const UniformVolume& volume, const FileFormatID format, const std::string& path )
{
  if ( ! volume.GetData() )
    {
    StdErr << "ERROR: cannot write volume that does not contain any data.\n";
    return;
    }

  DebugOutput( 3 ).GetStream().printf( "%s\nWriting %d x %d x %d voxels [%f x %f x %f mm total size].\n", path.c_str(), 
				       volume.GetDims()[0], volume.GetDims()[1], volume.GetDims()[2], volume.m_Size[0], volume.m_Size[1], volume.m_Size[2] );
  
  const TypedArray *data = volume.GetData();
  if ( data == NULL ) return;

  FileUtils::RecursiveMkPrefixDir( path );

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
    VolumeFromFile::WriteAnalyzeHdr( path, *actualVolume );
    break;
    }
    case FILEFORMAT_NIFTI_DETACHED: 
    case FILEFORMAT_NIFTI_SINGLEFILE: 
    {
    VolumeFromFile::WriteNifti( path, *actualVolume );
    break;
    }
    case FILEFORMAT_METAIMAGE: 
    {
    VolumeFromFile::WriteMetaImage( path, *actualVolume );
    break;
    }
    case FILEFORMAT_NRRD: 
    {
    VolumeFromFile::WriteNRRD( path, *actualVolume );
    break;
    }
    break;
    default:
      break;
    }
  
//  volume.SetMetaInfo( META_FS_PATH, path );
}

VolumeIO::Initializer::Initializer()
{
  if ( getenv( "IGS_WRITE_UNCOMPRESSED" ) || getenv( "CMTK_WRITE_UNCOMPRESSED" ) )
    VolumeIO::SetWriteCompressedOff();
}

VolumeIO::Initializer VolumeIO::Initializer::Instance;

} // namespace cmtk
