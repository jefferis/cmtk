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

#include "cmtkCompressedStream.h"

#include <bzlib.h>

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

CompressedStream::BZip2::BZip2( const char* filename ) 
{
  this->m_BzFile = BZ2_bzopen( filename, CMTK_FILE_MODE );
  if ( !this->m_BzFile ) 
    {
    throw 0;
    }
}

void 
CompressedStream::BZip2::Close()
{
  BZ2_bzclose( this->m_BzFile );
}

void
CompressedStream::BZip2::Rewind() 
{
  StdErr << "CompressedStream::BZip2::Rewind() is not implemented\n";
  throw ExitException( 1 );
}

size_t
CompressedStream::BZip2::Read ( void *data, size_t size, size_t count ) 
{
  const size_t bytesRead = BZ2_bzRead( &this->m_BzError, this->m_BzFile, data, size * count );
  this->m_BytesRead += bytesRead;
  return bytesRead / size;
}

bool
CompressedStream::BZip2::Get ( char &c)
{
  if ( this->Feof() || !this->Read( &c, sizeof(char), 1 ) )
    return false;

  return true;
}

int
CompressedStream::BZip2::Tell () const 
{
  return this->m_BytesRead;
}

bool
CompressedStream::BZip2::Feof () const 
{
  return this->m_BzError == BZ_STREAM_END;
}

} // namespace cmtk
