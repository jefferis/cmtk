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

#include "System/cmtkCompressedStream.h"

#include "System/cmtkConsole.h"
#include "System/cmtkMemory.h"

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>

#ifdef HAVE_SYS_TYPES_H
#  include <sys/types.h>
#endif

#ifdef _MSC_VER
#  include <process.h>
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
  : m_Reader( NULL )
{
  this->Open( filename );
}

CompressedStream::~CompressedStream () 
{
  this->Close();
}

bool
CompressedStream::Open ( const char *filename ) 
{
  this->Close();

  if ( Self::Stat( filename ) == 2 )
    {
    StdErr << "WARNING: file '" << filename << "' exists both compressed and uncompressed!\n";
    }
  
  const char *suffix = strrchr( filename, '.' );
  
  int compressed = 0;
  
  if ( suffix )
    for ( int i=0; ArchiveLookup[i].suffix && !compressed; ++i )
      compressed = compressed || ! strcmp( ArchiveLookup[i].suffix, suffix );

  try
    {
    if ( !compressed )
      {
      this->m_Reader = ReaderBase::SmartPtr( new Self::File( filename ) );
      }
    }
  catch (...)
    {
    }

  try 
    {
    if ( ! this->m_Reader )
      {
      bool result = false;
      for ( int i=0; ArchiveLookup[i].suffix && !result; ++i )
	result = this->OpenDecompressionPipe( filename, suffix, ArchiveLookup[i].command, ArchiveLookup[i].suffix );
      }
    }
  catch ( ... )
    {
    this->m_Reader = ReaderBase::SmartPtr( NULL );
    }
  
  return this->IsValid();
}

void
CompressedStream::Close()
{
  if ( this->m_Reader )
    {
    this->m_Reader->Close();
    this->m_Reader = ReaderBase::SmartPtr( NULL );
    }
}

bool
CompressedStream::OpenDecompressionPipe
( const char* filename, const char* suffix, const char* command, const char* compressedSuffix )
{
  char fname[PATH_MAX];

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
    if ( !strcmp( compressedSuffix, ".gz" ) ) 
      {
      this->m_Reader = ReaderBase::SmartPtr( new Self::Zlib( fname ) );
      }
#ifdef CMTK_USE_BZIP2
    else if ( !strcmp( compressedSuffix, ".bz2" ) ) 
      {
      this->m_Reader = ReaderBase::SmartPtr( new Self::BZip2( fname ) );
      }
#endif
    else
      {
      this->m_Reader = ReaderBase::SmartPtr( new Self::Pipe( filename, command ) );
      }
    }
  return this->IsValid();
}

std::string
CompressedStream::GetBaseName( const std::string& path )
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
