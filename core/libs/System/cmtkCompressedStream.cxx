/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
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

#include <cmtkCompressedStream.h>

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef _MSC_VER
#  include <process.h>
#endif

#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif

#include <cmtkConsole.h>
#include <cmtkMemory.h>

#if defined(_MSC_VER)
#define CMTK_FILE_MODE "rb"
#elif defined(__linux__)
#define CMTK_FILE_MODE "r"
#else
#define CMTK_FILE_MODE "rb"
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

const CompressedStream::ArchiveLookupEntry 
CompressedStream::ArchiveLookup[] = {
#ifdef _MSC_VER
  {".gz",  "gzip -cd %s > %s"},
  {".Z",   "gunzip -cd %s > %s"},
  {".bz",  "bzip -Q -cd %s > %s"},
  {".bz2", "bzip2 -cd %s > %s"},
#else
  {".gz",  "gzip -cd %s"},
  {".Z",   "gunzip -c %s"},
  {".bz",  "bzip -Q -cd %s"},
  {".bz2", "bzip2 -cd %s"},
#endif
  { NULL,   NULL} 
};

CompressedStream::CompressedStream ( const char *filename ) 
{
  File = NULL;
  GzFile = NULL;
  m_FilePointerMode = FILE_POINTER_INVALID;

  this->Open( filename );
}

CompressedStream::~CompressedStream () 
{
  this->Close();
}

bool CompressedStream::Open ( const char *filename ) 
{
  this->Close();
  BytesRead = 0;

  if ( Self::Stat( filename ) == 2 )
    {
    StdErr << "WARNING: file '" << filename << "' exists both compressed and uncompressed!\n";
    }
  
  const char *suffix = strrchr( filename, '.' );

  int compressed = 0;

  if ( suffix )
    for ( int i=0; ArchiveLookup[i].suffix && !compressed; ++i )
      compressed = compressed || ! strcmp( ArchiveLookup[i].suffix, suffix );

  if ( (!compressed) && ( File = fopen( filename, CMTK_FILE_MODE )) ) 
    { // file
    m_FilePointerMode = FILE_POINTER_FILE;
    } 
  else
    {
    bool result = false;
    for ( int i=0; ArchiveLookup[i].suffix && !result; ++i )
      result = this->OpenDecompressionPipe( filename, suffix, ArchiveLookup[i].command, ArchiveLookup[i].suffix );
    }
  return m_FilePointerMode != FILE_POINTER_INVALID;
}

void CompressedStream::Close()
{
  switch ( m_FilePointerMode ) 
    {
    case FILE_POINTER_FILE:
      fclose( File );
      File = NULL;
      break;
    case FILE_POINTER_PIPE:
#ifndef _MSC_VER
      pclose( File );
#else
      if ( m_FilePointerMode == FILE_POINTER_PIPE )
	remove( TempName );
#endif // # ifndef _MSC_VER
      File = NULL;
      break;
    case FILE_POINTER_ZLIB:
      gzclose( GzFile );
      GzFile = NULL;
      break;
    default:
      // no open stream - we're done.
      break;
    }

  m_FilePointerMode = FILE_POINTER_INVALID;
}

bool
CompressedStream::OpenDecompressionPipe
( const char* filename, const char* suffix, const char* command,
  const char* compressedSuffix )
{
  char fname[PATH_MAX], cmd[PATH_MAX];

  strcpy( fname, filename );
  if ( !suffix || strcmp( compressedSuffix, suffix ) )
    strcat( fname, compressedSuffix );

#ifdef _MSC_VER 
  for ( char *p=fname; *p; ++p )
    if ( *p == '/' ) *p = '\\';
#endif

  struct stat buf;
  if ( (! stat( fname, &buf )) && ( (buf.st_mode & S_IFREG) == S_IFREG ) ) 
    {
#ifndef _MSC_VER
    if ( !strcmp( compressedSuffix, ".gz" ) ) 
      {
      GzFile = gzopen( fname, CMTK_FILE_MODE );
      if ( GzFile ) 
	{
	m_FilePointerMode = FILE_POINTER_ZLIB;
	return true;
	}
      }
    
    if ( static_cast<size_t>( snprintf( cmd, sizeof( cmd ), command, fname ) ) >= sizeof( cmd ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in CompressedStream::OpenDecompressionPipe and will be truncated.\n";
      }
    errno = 0;

    File = popen( cmd, CMTK_FILE_MODE );
    if ( !File ) 
      {
      fprintf( stderr, "ERROR: popen() return NULL (errno=%d).\n", errno );
      perror( "System message" );
      return false;
      } 
    m_FilePointerMode = FILE_POINTER_PIPE;
    return true;
#else
    if ( snprintf( cmd, sizeof( cmd ), command, fname, tmpnam( TempName) ) >= sizeof( cmd ) )
      {
      StdErr << "WARNING: length of path exceeds system PATH_MAX in CompressedStream::OpenDecompressionPipe and will be truncated.\n";
      }
    
    _flushall();
    int sysReturn = system( cmd );

    if ( sysReturn ) 
      {
      fprintf( stderr, "Command %s returned %d\n", cmd, sysReturn );
      fprintf( stderr, "Errno = %d\n", errno );
      }
    
    File = fopen( TempName, CMTK_FILE_MODE);
    m_FilePointerMode = FILE_POINTER_FILE;
    return ( File != NULL );
#endif
  }
  return false;
}

int CompressedStream::Seek ( long int offset, int whence ) 
{
  switch ( m_FilePointerMode ) 
    {
    case FILE_POINTER_INVALID:
      return -1;
    case FILE_POINTER_FILE:
      return fseek( File, offset, whence );
    case FILE_POINTER_PIPE: 
    {
    char buffer[8192];
    int Result = 0;
    
    BytesRead += offset;
    while ( offset ) 
      if (offset<8192) 
	{
	Result += fread( buffer, sizeof(char), offset, File );
	offset=0;
	} 
      else
	{
	Result += fread( buffer, sizeof(char), 8192, File );
	offset-=8192;
	}
    return (Result == whence);
    }
    case FILE_POINTER_ZLIB:
      return gzseek( GzFile, offset, whence );
    default:
      return -1;
    }
  return -1;
}

size_t CompressedStream::Read ( void *data, size_t size, size_t count ) 
{
  BytesRead += size * count;
  if ( m_FilePointerMode == FILE_POINTER_ZLIB )
    return gzread( GzFile, data, size * count ) / size;
  return fread( data, size, count, File ); 
}

size_t
CompressedStream::Read ( void*& data ) 
{
  if ( m_FilePointerMode == FILE_POINTER_ZLIB ) 
    {
    gzseek( GzFile, 0, SEEK_END );
    size_t length = gztell( GzFile );
    gzseek( GzFile, 0, SEEK_SET );
    
    data = Memory::AllocateArray<char>( length );
    if ( length != static_cast<size_t>( gzread( GzFile, data, length ) ) ) 
      {
      Memory::Delete( data );
      data = NULL;
      length = 0;
      }
    
    return length;
    }

  fseek( File, 0, SEEK_END );
  size_t length = ftell( File );
  fseek( File, 0, SEEK_SET );
  
  data = Memory::AllocateArray<char>( length );
  if ( length != fread( data, 1, length, File ) ) 
    {
    Memory::Delete( data );
    data = NULL;
    length = 0;
    }
  
  return length;
}

int CompressedStream::Get ( char &c)
{
  int data;
  if ( m_FilePointerMode == FILE_POINTER_ZLIB )
    data = gzgetc( GzFile );
  else
    data = getc( File );

  if ( data != EOF ) 
    {
    ++BytesRead;
    c=(char) data;
    return 1;
    }
  return 0;
}

char*
CompressedStream::Gets ( char *const buffer, const int len )
{
  if ( m_FilePointerMode == FILE_POINTER_ZLIB )
    return gzgets( GzFile, buffer, len );
  else
    return fgets( buffer, len, File );
}

int CompressedStream::Tell () const 
{
  switch ( m_FilePointerMode ) 
    {
    case FILE_POINTER_FILE:
      return ftell( File );
    case FILE_POINTER_PIPE:
      return BytesRead;
    case FILE_POINTER_ZLIB:
      return gztell( GzFile );
    default:
      return -1;
    }
  return -1;
}

int CompressedStream::Feof () const 
{
  if ( m_FilePointerMode == FILE_POINTER_ZLIB )
    return gzeof( GzFile );
  return feof( File );
}

std::string CompressedStream::GetBaseName( const std::string& path )
{
  const size_t suffixPos = path.rfind( '.' );
  
  if ( suffixPos != std::string::npos ) 
    {
    for ( int i = 0; ArchiveLookup[i].suffix; ++i )
      {
      const size_t suffixLen = strlen( ArchiveLookup[i].suffix );
      if ( !path.compare( suffixPos, suffixLen, ArchiveLookup[i].suffix, suffixLen ) )
	{
	return path.substr( 0, suffixPos );
	}
      }
    }
  return path;
}

int 
CompressedStream::Stat( const char *path, struct stat* buf )
{
  std::string baseName = CompressedStream::GetBaseName( path );

  struct stat statbuf;
  if ( ! buf )
    buf = &statbuf;

  const bool existsUncompressed = ! stat( baseName.c_str(), buf );
  
  for ( int i = 0; ArchiveLookup[i].suffix; ++i ) 
    {
    const std::string cpath = baseName + std::string( ArchiveLookup[i].suffix );
    if ( ! stat( cpath.c_str(), buf ) ) 
      return existsUncompressed ? 2 : 1;
    }
  
  return existsUncompressed ? 0 : -1;
}

} // namespace cmtk
