/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2014 SRI International
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

#include <System/cmtkCompressedStream.h>

#include <System/cmtkConsole.h>
#include <System/cmtkMemory.h>
#include <System/cmtkMountPoints.h>

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
  {".Z",    "gunzip -cd %s > %s"},
  {".gz",   "gzip -cd %s > %s"},
  {".bz",   "bzip -Q -cd %s > %s"},
  {".bz2",  "bzip2 -cd %s > %s"},
  {".lzma", "xz -cd %s > %s"},
  {".xz",   "xz -cd %s > %s"},
#else
  {".Z",    "gunzip -c %s"},
  {".gz",   "gzip -cd %s"},
  {".bz",   "bzip -Q -cd %s"},
  {".bz2",  "bzip2 -cd %s"},
  {".lzma", "xz -cd %s"},
  {".xz",   "xz -cd %s"},
#endif
  { NULL,   NULL} 
};

CompressedStream::CompressedStream ( const std::string& filename ) 
  : m_Reader( NULL ),
    m_Compressed( false )
{
  this->Open( MountPoints::Translate( filename ) );
}

CompressedStream::~CompressedStream () 
{
  this->Close();
}

bool
CompressedStream::Open ( const std::string& filename ) 
{
  this->Close();

  if ( Self::Stat( filename.c_str() ) == 2 )
    {
    StdErr << "WARNING: file '" << filename << "' exists both compressed and uncompressed!\n";
    }
  
  this->m_Compressed = false;
  
  std::string suffix = "";
  const size_t period = filename.rfind( '.' );
  if ( period != std::string::npos )
    {
    suffix = filename.substr( period, std::string::npos );
    for ( int i=0; ArchiveLookup[i].suffix && !this->m_Compressed; ++i )
      this->m_Compressed = this->m_Compressed || ( suffix == ArchiveLookup[i].suffix );
    }
  
  try
    {
    if ( !this->m_Compressed )
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
      this->m_Compressed = true;
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
( const std::string& filename, const std::string& suffix, const char* command, const char* compressedSuffix )
{
  std::string fname = filename;
  if ( suffix != compressedSuffix )
    fname = fname + compressedSuffix;

#ifdef _MSC_VER 
  std::replace( fname.begin(), fname.end(), '/', '\\' );
#endif

  Self::StatType buf;
  if ( (! 
#ifdef CMTK_USE_STAT64
          stat64( fname.c_str(), &buf )
#else
          stat( fname.c_str(), &buf )
#endif
       ) && ( (buf.st_mode & S_IFREG) == S_IFREG ) ) 
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
#ifdef CMTK_USE_LZMA
    else if ( !strcmp( compressedSuffix, ".lzma" ) ) 
      {
      this->m_Reader = ReaderBase::SmartPtr( new Self::LZMA( fname ) );
      }
#endif
    else
      {
      this->m_Reader = ReaderBase::SmartPtr( new Self::Pipe( fname, command ) );
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
    const std::string fileSuffix = path.substr( suffixPos, std::string::npos );
    for ( int i = 0; ArchiveLookup[i].suffix; ++i )
      {
      if ( fileSuffix == ArchiveLookup[i].suffix )
	{
	return path.substr( 0, suffixPos );
	}
      }
    }
  return path;
}

int 
CompressedStream::Stat( const std::string& path, Self::StatType* buf )
{
  const std::string baseName = CompressedStream::GetBaseName( MountPoints::Translate( path ) );

  Self::StatType statbuf;
  if ( ! buf )
    buf = &statbuf;

#ifdef CMTK_USE_STAT64
  const bool existsUncompressed = ! stat64( baseName.c_str(), buf );
#else
  const bool existsUncompressed = ! stat( baseName.c_str(), buf );
#endif
  
  for ( int i = 0; ArchiveLookup[i].suffix; ++i ) 
    {
    const std::string cpath = baseName + std::string( ArchiveLookup[i].suffix );
    if ( ! 
#ifdef CMTK_USE_STAT64
           stat64( cpath.c_str(), buf )
#else
           stat( cpath.c_str(), buf )
#endif
        ) 
      return existsUncompressed ? 2 : 1;
    }
  
  return existsUncompressed ? 0 : -1;
}

} // namespace cmtk
