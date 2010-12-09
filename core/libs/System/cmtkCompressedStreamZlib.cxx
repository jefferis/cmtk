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

#include <zlib.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

CompressedStream::Zlib::Zlib( const char* filename )
{
  this->m_GzFile = gzopen( filename, CMTK_FILE_MODE );
  if ( !this->m_GzFile ) 
    {
    throw 0;
    }
}

void 
CompressedStream::Zlib::Close()
{
  gzclose( this->m_GzFile );
}

void
CompressedStream::Zlib::Rewind ()
{
  gzrewind( this->m_GzFile );
  this->CompressedStream::ReaderBase::Rewind();
}

int
CompressedStream::Zlib::Seek ( const long int offset, int whence ) 
{
  return gzseek( this->m_GzFile, offset, whence );
}

size_t
CompressedStream::Zlib::Read ( void *data, size_t size, size_t count ) 
{
  return gzread( this->m_GzFile, data, size * count ) / size;
}

bool
CompressedStream::Zlib::Get ( char &c)
{
  const int data = gzgetc( this->m_GzFile );
  if ( data != EOF ) 
    {
    c=(char) data;
    return true;
    }

  return false;
}

int
CompressedStream::Zlib::Tell () const 
{
  return gztell( this->m_GzFile );
}

bool
CompressedStream::Zlib::Feof () const 
{
  return gzeof( this->m_GzFile );
}

} // namespace cmtk
