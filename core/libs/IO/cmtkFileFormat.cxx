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

#include <cmtkFileFormat.h>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#include <stdio.h>
#include <string.h>
#include <limits.h>

#include <cmtkTypes.h>
#include <cmtkCompressedStream.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/// Structure that holds information on magic file number.
typedef 
struct
{
  /// Offset of the magic number within the file.
  unsigned short offset;
  /// Magic number as string of characters.
  const char *magicString;
  /// Length of magic number string (necessary due to Null characters).
  const size_t magicStringLength;
} FileFormatMagic;

/// Magic number records for known file types.
const FileFormatMagic FileFormatMagicNumbers[] = {
  { 0, NULL, 0 }, // NEXIST
  { 0, NULL, 0 }, // STUDY
  { 0, NULL, 0 }, // STUDYLIST
  { 0, "! TYPEDSTREAM", 13 }, // TypedStream archive
  { 0, "P5", 2 }, // PGM
  { 128, "DICM", 4 }, // DICOM
  { 0, "Modality :=", 11 }, // VANDERBILT
  { 252, "\x2A\xE3\x89\xB8", 4 }, // GIPL
  { 0, "# AmiraMesh 3D", 14 },
  { 0, NULL, 0 }, // RAW
  { 0, NULL, 0 }, // RAW3D
  { 0, "W\\YX", 4 }, // ACCURAY
  { 54, "\x39\x30", 2 }, // BIORAD
  { 344, "ni1\x00", 4 }, // Nifti, detached header
  { 344, "n+1\x00", 4 }, // Nifti, single file
  { 0, "AVW_ImageFile", 13 }, // Analyze AVW file.
  { 0, NULL, 0 }, // MetaImage
  { 0, "NRRD", 4 }, //NRRD
  { 0, "\x5C\x01\x00\x00", 4 }, // Analyze little endian
  { 0, "\x00\x00\x01\x5C", 4 }, // Analyze big endian
  { 0, NULL, 0 }  // Unknown.
};

const char* FileFormatName[] =
{
  /// File of this name does not exist.
  "",
  /// Path is a typedstream study directory.
  "STUDY",
  /// Path is a typedstream studylist directory.
  "STUDYLIST",
  /// Path is a typedstream archive.
  "TYPEDSTREAM",
  /// Path is a PGM image.
  "PGM",
  /// Path is a DICOM image.
  "DICOM",
  /// Path is a Vanderbilt image description file.
  "VANDERBILT",
  /// Path is a GIPL (Guy's Hospital Image Processing Lab) image file.
  "GIPL",
  /// Path is a Amira image file.
  "AMIRA",
  /// Path is some raw binary file (2-D).
  "RAW-DATA",
  /// Path is some raw binary file (3-D).
  "RAW3D",
  /// Path is an Accuray file.
  "ACCURAY",
  /// Path is a BioRad .PIC file.
  "BIORAD",
  /// Nifti, detached header
  "NIFTI-DETACHED-HEADER",
  /// Nifti, single file
  "NIFTI-SINGLE-FILE",
  /// Path is an Analyze AVW file.
  "ANALYZE-AVW",
  /// MetaImage
  "METAIMAGE",
  /// NRRD
  "NRRD",
  /// Path is an Analyze 7.5 file in little endian.
  "ANALYZE-HDR-LITTLEENDIAN",
  /// Path is an Analyze 7.5 file in big endian.
  "ANALYZE-HDR-BIGENDIAN",
  /** File type cannot be determined.
   * This ID always has to be the last one!
   */
  NULL
};

FileFormatID 
FileFormat::Identify( const char* path )
{
  struct stat buf;
  if ( CompressedStream::Stat( path, &buf ) < 0 ) 
    return FILEFORMAT_NEXIST;

  if ( buf.st_mode & S_IFDIR ) 
    return FileFormat::IdentifyDirectory( path );
  else if ( buf.st_mode & S_IFREG ) 
    return FileFormat::IdentifyFile( path );

  return FILEFORMAT_NEXIST;
}

FileFormatID 
FileFormat::GetID( const char* name )
{
  if ( name ) 
    {
    for ( unsigned int idx = 0; FileFormatName[idx]; ++idx ) 
      {
      if ( ! strcmp( FileFormatName[idx], name ) )
	return static_cast<FileFormatID>( idx );
      }
    }
  return FILEFORMAT_UNKNOWN;
}

const char* 
FileFormat::Describe( const FileFormatID id )
{
  switch ( id ) 
    {
    case FILEFORMAT_NEXIST:
      return "File or directory does not exist.";
    case FILEFORMAT_STUDY:
      return "Typedstream study archive [Directory].";
    case FILEFORMAT_STUDYLIST:
      return "Typedstream studylist archive [Directory].";
    case FILEFORMAT_PGM:
      return "PGM image file [File].";
    case FILEFORMAT_DICOM:
      return "DICOM image file [File].";
    case FILEFORMAT_VANDERBILT:
      return "Vanderbilt header/image file combination [File].";
    case FILEFORMAT_GIPL:
      return "GIPL image file [File].";
    case FILEFORMAT_AMIRA:
      return "AmiraMesh image file [File].";
    case FILEFORMAT_ACCURAY:
      return "Accuray image file [File].";
    case FILEFORMAT_BIORAD:
      return "BioRad image file [File].";
    case FILEFORMAT_NIFTI_DETACHED:
      return "NIFTI detached header+image [File]";      
    case FILEFORMAT_NIFTI_SINGLEFILE:
      return "NIFTI single file [File]";      
    case FILEFORMAT_ANALYZE_HDR:
      return "Analyze 7.5 file [Header+Binary File/Little Endian].";
    case FILEFORMAT_ANALYZE_HDR_BIGENDIAN:
      return "Analyze 7.5 file [Header+Binary File/Big Endian].";
    case FILEFORMAT_ANALYZE_AVW:
      return "Analyze AVW file [File].";
    case FILEFORMAT_RAW:
      return "RAW image file [File].";
    case FILEFORMAT_NRRD:
      return "Nrrd image file [File].";
    case FILEFORMAT_UNKNOWN:
    default:
      return "Unknown format.";
    }
  return "ILLEGAL ID tag in FileFormat::Describe().";
}

FileFormatID 
FileFormat::IdentifyDirectory( const char* path )
{
  char filename[PATH_MAX];
  struct stat buf;

  snprintf( filename, sizeof( filename ), "%s/images", path );
  if ( (!stat( filename, &buf )) && ( buf.st_mode & S_IFREG ) )
    return FILEFORMAT_STUDY;

  snprintf( filename, sizeof( filename ), "%s/images.gz", path );
  if ( (!stat( filename, &buf )) && ( buf.st_mode & S_IFREG ) )
    return FILEFORMAT_STUDY;

  snprintf( filename, sizeof( filename ), "%s/studylist", path );
  if ( (!stat( filename, &buf )) && ( buf.st_mode & S_IFREG ) )
    return FILEFORMAT_STUDYLIST;

  snprintf( filename, sizeof( filename ), "%s/studylist.gz", path );
  if ( (!stat( filename, &buf )) && ( buf.st_mode & S_IFREG ) )
    return FILEFORMAT_STUDYLIST;

  return FILEFORMAT_UNKNOWN;
}

FileFormatID 
FileFormat::IdentifyFile( const char* path )
{
  CompressedStream stream( path );
  if ( ! stream.IsValid() )
    return FILEFORMAT_NEXIST;
    
  char buffer[348];
  memset( buffer, 0, sizeof( buffer ) );
  stream.Read( buffer, 1, 348 );
  
  FileFormatID id = FILEFORMAT_NEXIST; 
  while ( id != FILEFORMAT_UNKNOWN ) 
    {
    if ( FileFormatMagicNumbers[id].magicString ) 
      {
      if ( !memcmp( buffer+FileFormatMagicNumbers[id].offset, FileFormatMagicNumbers[id].magicString, FileFormatMagicNumbers[id].magicStringLength ) )
	return id;
      }
    char cid = static_cast<char>( id );
    id = static_cast<FileFormatID>( ++cid );
    }
  
  return FILEFORMAT_UNKNOWN;
}

} // namespace cmtk
