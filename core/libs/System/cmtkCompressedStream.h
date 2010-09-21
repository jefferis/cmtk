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

#ifndef __cmtkCompressedStream_h_included_
#define __cmtkCompressedStream_h_included_

#include <cmtkconfig.h>

#include "System/cmtkCannotBeCopied.h"

#include "System/cmtkException.h"

#include <stdio.h>

#include <zlib.h>

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#include <string>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Stream with on-the-fly decompression
class CompressedStream :
  /// Make class uncopyable via inheritance.
  private CannotBeCopied
{
public:
  /// This class.
  typedef CompressedStream Self;
  
  /// Create stream object without opening any files.
  CompressedStream() : m_FilePointerMode( FILE_POINTER_INVALID ) {}

  /// Create stream from filename.
  CompressedStream ( const char *filename );
  
  /// Dispose stream object.
  ~CompressedStream ();

  /// Return validity of stream.
  bool IsValid() const { return m_FilePointerMode != FILE_POINTER_INVALID; }

  /// Open new stream from filename.
  bool Open( const char *filename );
  
  /// Close current file stream.
  void Close();

  /** Set filepointer.
   * If the object represents a pipe with on-the-fly decompression, only
   * the set mode "SEEK_CUR" is supported and will be simulated by successive
   * reading of 8kB data blocks from the pipe.
   *@param offset Offset the file pointer is set to, depending on the value of
   * whence.
   *@param whence File pointer set mode as defined for fseek.
   */
  int Seek ( long int offset, int whence );

  /// Read block of data.
  size_t Read ( void *data, size_t size, size_t count );

  /// Read a single character from the stream.
  int Get ( char &c);

  /// Get string function
  char *Gets ( char *const buffer, const int len );

  /// Return number of bytes read from stream.
  int Tell () const;

  /// Return 1 if and only if end of file reached.
  int Feof () const;

  /** Return base name of a path without compression suffix.
   */
  static std::string GetBaseName( const std::string& path );

  /** Do stat() on compressed file.
   *\return -1 if the file does not exist; 0 if the file exists under the given name; 1 if the file exists
   * with an additional suffix corresponding to one of the supported compression schemes, e.g., ".gz"; 2 if
   * the file exists with both its plain name AND at least one compressed suffix. The last case indicates a
   * potential consistency problem because it is not clear, which file should be read.
   */
  static int Stat( const char *path, struct stat *const buf = NULL );

private:
  /// File pointer to read from.
  FILE *File;

  /// GZIP file pointer when using zlib decompression.
  gzFile GzFile;

  /// Enum for different internal file modes.
  typedef enum 
  {
    /// All file pointers are invalid.
    FILE_POINTER_INVALID,
    /// File pointer points to open regular file.
    FILE_POINTER_FILE,
    /// File pointer points to open pipe.
    FILE_POINTER_PIPE,
    /// GzFile pointer points to open zlib file object.
    FILE_POINTER_ZLIB
  } FilePointerMode;

  /// Flag indicating whether stream is file, pipe, or zlib file.
  FilePointerMode m_FilePointerMode;

#ifdef _MSC_VER
  /** Temporary filename.
   * On platforms where no on-the-fly decompression by via a pipe is possible
   * (i.e. Windows), this holds the name of a temporary file used to hold
   * the decompressed input file.
   */
  char TempName[256];
#endif

  /// Count number of bytes read from file or pipe.
  long int BytesRead;

  /** Open decompressing pipe.
   * A suffix is appended to the desired filename, unless the name has
   * has already been given this suffix. Afterwards, a pipe through a
   * user-defined command is opened.
   *@param filename The name of the file to read through the pipe.
   *@param command The command to pipe the file through. This should be the
   * name of a program that reads the file with the given name (referenced by
   * "%s" in the command string) and writes its output to stdout.
   *@param suffix The file suffix corresponding to input files for the given
   * pipe command, e.g. ".gz" for gzip.
   */
  bool OpenDecompressionPipe ( const char* filename, const char* suffix, const char* command, const char* compressedSuffix );

  /** Entry for the suffix-to-archiver assignment table.
   */
  typedef struct 
  {
    /// A compressed file suffix such as .gz
    const char *suffix;
    
    /** Command to decompress a file with this suffix to a pipe.
     * For .gz, this command would be "gzip -df".
     */
    const char *command;
  } ArchiveLookupEntry;
  
  /// The suffix-to-archiver assignment table.
  static const ArchiveLookupEntry ArchiveLookup[];
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkCompressedStream_h_included_
