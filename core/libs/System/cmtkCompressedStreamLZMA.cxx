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

#include <lzmadec.h>

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

CompressedStream::LZMA::LZMA( const char* filename )
{
  this->m_File = lzmadec_open( filename );
  if ( !this->m_File ) 
    {
    StdErr << "ERROR: lzmadec_open() failed for file '" << filename << "'\n";
    throw ExitException( 1 );
    }
}

void 
CompressedStream::LZMA::Close()
{
  lzmadec_close( this->m_File );
}

size_t
CompressedStream::LZMA::Read ( void *data, size_t size, size_t count ) 
{
  const size_t result = lzmadec_read( this->m_File, reinterpret_cast<uint8_t*>( data ), size * count );
  this->m_BytesRead += result;
  return result / size;  
}

bool
CompressedStream::LZMA::Get ( char &c)
{
  const int data = lzmadec_getc( this->m_File );
  if ( data != EOF ) 
    {
    c=(char) data;
    ++this->m_BytesRead;
    return true;
    }

  return false;
}

void
CompressedStream::LZMA::Rewind ()
{
  lzmadec_rewind( this->m_File );
  this->CompressedStream::ReaderBase::Rewind();
}

int
CompressedStream::LZMA::Tell () const 
{
  return lzmadec_tell( this->m_File );
}

bool
CompressedStream::LZMA::Feof () const 
{
  return lzmadec_eof( this->m_File );
}

} // namespace cmtk
