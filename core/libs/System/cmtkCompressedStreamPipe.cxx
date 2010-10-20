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

#include <System/cmtkCompressedStream.h>

#include <System/cmtkConsole.h>

#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <errno.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

CompressedStream::Pipe::Pipe( const char* filename, const char* command )
{
  char cmd[PATH_MAX];

#ifndef _MSC_VER  
  if ( static_cast<size_t>( snprintf( cmd, sizeof( cmd ), command, filename ) ) >= sizeof( cmd ) )
    {
    StdErr << "WARNING: length of path exceeds system PATH_MAX in CompressedStream::OpenDecompressionPipe and will be truncated.\n";
    }
  errno = 0;
  
  this->m_File = popen( cmd, CMTK_FILE_MODE );
  if ( !this->m_File ) 
    {
    fprintf( stderr, "ERROR: popen() return NULL (errno=%d).\n", errno );
    perror( "System message" );
    throw 0;
    } 
#else
  if ( snprintf( cmd, sizeof( cmd ), command, filename, tmpnam( this->m_TempName) ) >= sizeof( cmd ) )
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
  
  this->m_File = fopen( this->m_TempName, CMTK_FILE_MODE);
  
  if ( !this->m_File )
    {
    throw 0;
    }
#endif

  this->m_BytesRead = 0;
}

void 
CompressedStream::Pipe::Close()
{
#ifndef _MSC_VER
  pclose( this->m_File );
#else
  remove( this->m_TempName );
#endif // # ifndef _MSC_VER
}

int
CompressedStream::Pipe::Seek ( long int offset, int whence ) 
{
  char buffer[Self::BlockSize];
  int result = 0;
  
  this->m_BytesRead += offset;
  while ( offset > 0 ) 
    {
    if ( static_cast<size_t>( offset ) < Self::BlockSize ) 
      {
      result += fread( buffer, sizeof(char), offset, this->m_File );
      offset=0;
      } 
    else
      {
      result += fread( buffer, sizeof(char), Self::BlockSize, this->m_File );
      offset -= Self::BlockSize;
      }
    }
  return (result == whence);
}

size_t
CompressedStream::Pipe::Read( void *data, size_t size, size_t count ) 
{
  return fread( data, size, count, this->m_File );
}

bool
CompressedStream::Pipe::Get ( char &c)
{
  const int data = fgetc( this->m_File );
  if ( data != EOF ) 
    {
    c=(char) data;
    return true;
    }

  return false;
}

int
CompressedStream::Pipe::Tell () const 
{
  return this->m_BytesRead;
}

bool
CompressedStream::Pipe::Feof () const 
{
  return feof( this->m_File );
}

} // namespace cmtk
