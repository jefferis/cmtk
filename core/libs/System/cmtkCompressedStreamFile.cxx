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

#include <stdlib.h>
#include <stdio.h>
#include <errno.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

CompressedStream::File::File( const char* filename )
{
  this->m_File = fopen( filename, CMTK_FILE_MODE );
  if ( !this->m_File ) 
    {
    throw 0;
    }
}

void 
CompressedStream::File::Close()
{
  fclose( this->m_File );
}

int
CompressedStream::File::Seek ( long int offset, int whence ) 
{
  return fseek( this->m_File, offset, whence );
}

size_t
CompressedStream::File::Read( void *data, size_t size, size_t count ) 
{
  return fread( data, size, count, this->m_File );
}

bool
CompressedStream::File::Get ( char &c)
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
CompressedStream::File::Tell () const 
{
  return ftell( this->m_File );
}

bool
CompressedStream::File::Feof () const 
{
  return feof( this->m_File );
}

} // namespace cmtk
