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

#include <bzlib.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

CompressedStreamBZip2BZip2( const char* filename )
{
  this->m_BzFile = BZ2_bzopen( filename, CMTK_FILE_MODE );
  if ( !this->m_BzFile ) 
    {
    throw 0;
    }
}

void 
CompressedStreamBZip2::Close()
{
  BZ2_bzclose( this->m_BzFile );
}

int
CompressedStreamBZip2::Seek ( long int offset, int whence ) 
{
  return gzseek( this->m_BzFile, offset, whence );
}

size_t
CompressedStreamBZip2::Read ( void *data, size_t size, size_t count ) 
{
  const size_t bytesRead = BZ2_bzread( this->m_BzFile, data, size * count );
  this->m_BytesRead += bytesRead;
  return bytesRead / size;
}

bool
CompressedStreamBZip2::Get ( char &c)
{
  return false;
}

char*
CompressedStreamBZip2::Gets ( char *const buffer, const int len )
{
  return NULL;
}

int
CompressedStreamBZip2::Tell () const 
{
  return this->m_BytesRead;
}

bool
CompressedStreamBZip2::Feof () const 
{
  return bzerror == BZ_STREAM_END;
}

} // namespace cmtk
