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

#ifndef __cmtkCompressedStream_h_included_
#define __cmtkCompressedStream_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>

#include <System/cmtkException.h>
#include <System/cmtkSmartPtr.h>

#include <stdio.h>

#include <zlib.h>

#ifdef CMTK_USE_BZIP2
#  include <bzlib.h>
#endif

#ifdef CMTK_USE_LZMA
#  include <lzmadec.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

#include <string>

#if defined(_MSC_VER)
#define CMTK_FILE_MODE "rb"
#else
#define CMTK_FILE_MODE "r"
#endif

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
  CompressedStream() : m_Reader( NULL ) {};
  
  /// Create stream from filename.
  CompressedStream ( const char *filename );
  
  /// Dispose stream object.
  ~CompressedStream ();
  
  /// Return validity of stream.
  bool IsValid() const
  {
    return (this->m_Reader != NULL); 
  }
  
  /// Return flag whether this stream is actually using decompression or simply reads a file as-is.
  bool IsCompressed() const
  {
    return this->m_Compressed;
  }
  
  /// Open new stream from filename.
  bool Open( const char *filename );
  
  /// Close current file stream.
  void Close();
  
  /** Set filepointer to beginning of file.
   */
  void Rewind()
  {
    this->m_Reader->Rewind();
  }

  /** Set filepointer.
   * If the object represents a pipe with on-the-fly decompression, only
   * the set mode "SEEK_CUR" is supported and will be simulated by successive
   * reading of 8kB data blocks from the pipe.
   *\param offset Offset the file pointer is set to, depending on the value of
   * whence.
   *\param whence File pointer set mode as defined for fseek.
   */
  int Seek ( const long int offset, int whence )
  {
    return this->m_Reader->Seek( offset, whence );
  }

  /// Read block of data.
  size_t Read ( void *data, size_t size, size_t count )
  {
    return this->m_Reader->Read( data, size, count );
  }

  /// Read a single character from the stream.
  int Get ( char &c)
  {
    return this->m_Reader->Get( c );
  }

  /// Return number of bytes read from stream.
  int Tell () const
  {
    return this->m_Reader->Tell();
  }

  /// Return 1 if and only if end of file reached.
  bool Feof () const
  {
    return this->m_Reader->Feof();
  }

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
  /** Open decompressing pipe.
   * A suffix is appended to the desired filename, unless the name has
   * has already been given this suffix. Afterwards, a pipe through a
   * user-defined command is opened.
   *\param filename The name of the file to read through the pipe.
   *\param suffix Actual suffix of the file to be read through the pipe.
   *\param command The command to pipe the file through. This should be the
   * name of a program that reads the file with the given name (referenced by
   * "%s" in the command string) and writes its output to stdout.
   *\param compressedSuffix The file suffix corresponding to input files for the given
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

private:
  /// Abstract base class for low-level reader engines.
  class ReaderBase
  {
  public:
    /// This class.
    typedef ReaderBase Self;
    
    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Default constructor.
    ReaderBase() : m_BytesRead( 0 ) {}

    /// Virtual destructor.
    virtual ~ReaderBase() {}

    /// Close current file stream.
    virtual void Close() = 0;

    /// Reset read pointer to beginning of stream.
    virtual void Rewind()
    {
      this->m_BytesRead = 0;
    }
    
    /** Set filepointer.
      * This class implements a naive seek that optionally calls "this->Rewind()" (if
      * "whence" is SEEK_SET, then reads "offset" bytes from the input to position
      * read pointer.
      *\param offset Offset the file pointer is set to, depending on the value of
      * whence.
      *\param whence File pointer set mode as defined for fseek.
      */
    virtual int Seek ( const long int offset, int whence );
    
    /// Read block of data.
    virtual size_t Read ( void *data, size_t size, size_t count ) = 0;
    
    /// Read a single character from the stream.
    virtual bool Get ( char &c) = 0;
    
    /// Return number of bytes read from stream.
    virtual int Tell () const = 0;
    
    /// Return 1 if and only if end of file reached.
    virtual bool Feof () const = 0;

  private:
    /// Block size for fake seek() operation.
    static const size_t SeekBlockSize = 8192;

  protected:
    /// Count number of bytes read from file or pipe.
    size_t m_BytesRead;
  };

  /// Class for uncompressed file-based reader engine.
  class File
    : public ReaderBase
  {
  public:
    /// This class.
    typedef File Self;
    
    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Open new stream from filename.
    File( const char *filename );
    
    /// Virtual destructor.
    virtual ~File() {}

    /// Close current file stream.
    virtual void Close();
    
    /// Reset read pointer to beginning of stream.
    virtual void Rewind();
    
    /** Set filepointer.
      * If the object represents a pipe with on-the-fly decompression, only
      * the set mode "SEEK_CUR" is supported and will be simulated by successive
      * reading of 8kB data blocks from the pipe.
      *\param offset Offset the file pointer is set to, depending on the value of
      * whence.
      *\param whence File pointer set mode as defined for fseek.
      */
    virtual int Seek ( const long int offset, int whence );
    
    /// Read block of data.
    virtual size_t Read ( void *data, size_t size, size_t count );
    
    /// Read a single character from the stream.
    virtual bool Get ( char &c);
    
    /// Return number of bytes read from stream.
    virtual int Tell () const;
    
    /// Return 1 if and only if end of file reached.
    virtual bool Feof () const;

  private:
    /// File pointer.
    FILE* m_File;    
  };

  /// Class for reader engine using pipe.
  class Pipe
    : public ReaderBase
  {
  public:
    /// This class.
    typedef Pipe Self;
    
    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Open new pipe from filename.
    Pipe( const char* filename, const char* command );
    
    /// Virtual destructor.
    virtual ~Pipe() {}

    /// Close current file stream.
    virtual void Close();
    
    /// Reset read pointer to beginning of stream.
    virtual void Rewind();
    
    /// Read block of data.
    virtual size_t Read ( void *data, size_t size, size_t count );
    
    /// Read a single character from the stream.
    virtual bool Get ( char &c);
    
    /// Return number of bytes read from stream.
    virtual int Tell () const;
    
    /// Return 1 if and only if end of file reached.
    virtual bool Feof () const;

  private:
    /// File pointer.
    FILE* m_File;    
    
#ifdef _MSC_VER
    /** Temporary filename.
   * On platforms where no on-the-fly decompression by via a pipe is possible
   * (i.e. Windows), this holds the name of a temporary file used to hold
   * the decompressed input file.
   */
    char m_TempName[256];
#endif
  };

#ifdef CMTK_USE_BZIP2
  /// Class for BZip2-based reader engine.
  class BZip2 
    : public ReaderBase
  {
  public:
    /// This class.
    typedef BZip2 Self;
    
    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Open new stream from filename.
    BZip2( const char *filename );
    
    /// Virtual destructor.
    virtual ~BZip2() {}

    /// Close current file stream.
    virtual void Close();
    
    /// Reset read pointer to beginning of stream.
    virtual void Rewind();
    
    /// Read block of data.
    virtual size_t Read ( void *data, size_t size, size_t count );
    
    /// Read a single character from the stream.
    virtual bool Get ( char &c);
    
    /// Return number of bytes read from stream.
    virtual int Tell () const;
    
    /// Return 1 if and only if end of file reached.
    virtual bool Feof () const;

  private:
    /// BZip2 file pointer.
    BZFILE* m_BzFile;    

    /// BZip2 error variable.
    int m_BzError;
  };
#endif

#ifdef CMTK_USE_LZMA
  /// Class for Zlib-based reader engine.
  class LZMA 
    : public ReaderBase
  {
  public:
    /// This class.
    typedef LZMA Self;
    
    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Open new stream from filename.
    LZMA( const char *filename );
    
    /// Virtual destructor.
    virtual ~LZMA() {}

    /// Close current file stream.
    virtual void Close();
    
    /// Reset read pointer to beginning of stream.
    virtual void Rewind();
    
    /// Read block of data.
    virtual size_t Read ( void *data, size_t size, size_t count );
    
    /// Read a single character from the stream.
    virtual bool Get ( char &c);
    
    /// Return number of bytes read from stream.
    virtual int Tell () const;
    
    /// Return 1 if and only if end of file reached.
    virtual bool Feof () const;

  private:
    /// LZMA file pointer.
    lzmadec_FILE* m_File;    
  };
#endif // #ifdef CMTK_USE_LZMA

  /// The low-level reader object.
  ReaderBase::SmartPtr m_Reader;

  /// Flag whether current stream is from a compressed source.
  bool m_Compressed;

public:
  /// Class for Zlib-based reader engine.
  class Zlib 
    : public ReaderBase
  {
  public:
    /// This class.
    typedef Zlib Self;
    
    /// Smart pointer to this class.
    typedef SmartPointer<Self> SmartPtr;

    /// Open new stream from filename.
    Zlib( const char *filename );
    
    /// Virtual destructor.
    virtual ~Zlib() {}

    /// Close current file stream.
    virtual void Close();
    
    /// Reset read pointer to beginning of stream.
    virtual void Rewind();
    
    /** Set filepointer.
      * If the object represents a pipe with on-the-fly decompression, only
      * the set mode "SEEK_CUR" is supported and will be simulated by successive
      * reading of 8kB data blocks from the pipe.
      *\param offset Offset the file pointer is set to, depending on the value of
      * whence.
      *\param whence File pointer set mode as defined for fseek.
      */
    virtual int Seek ( const long int offset, int whence );
    
    /** Read block of data.
     * This function, unlike gzread(), is safe to use on very large files.
     *\see http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=587912
     *\remark Yaroslav Halchenko <debian@onerussian.com> pointed this out.
     */

    virtual size_t Read ( void *data, size_t size, size_t count );
    
    /** Safe, static write function.
     * This is simply here as a service to other classes, which may use this as a plug-in replacement
     * for gzwrite() to work around a large-file problem in zlib.
     *\see http://bugs.debian.org/cgi-bin/bugreport.cgi?bug=587912
     *\remark Yaroslav Halchenko <debian@onerussian.com> pointed this out.
     */
    static size_t StaticSafeWrite ( gzFile file, const void *data, size_t size );

    /// Read a single character from the stream.
    virtual bool Get ( char &c );
    
    /// Return number of bytes read from stream.
    virtual int Tell () const;
    
    /// Return 1 if and only if end of file reached.
    virtual bool Feof () const;

  private:
    /// Zlib file pointer.
    gzFile m_GzFile;    
  };
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkCompressedStream_h_included_
