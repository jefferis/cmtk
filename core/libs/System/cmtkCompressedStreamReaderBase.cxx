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

int
cmtk::CompressedStream::ReaderBase::Seek( const long int offset, int whence ) 
{
  if ( whence == SEEK_SET )
    this->Rewind();
  
  char buffer[Self::SeekBlockSize];
  int result = 0;
  
  for ( int stillToRead = offset; stillToRead > 0; ) 
    {
    if ( static_cast<size_t>( stillToRead ) < Self::SeekBlockSize )
      {
      result += this->Read( buffer, sizeof(char), stillToRead );
      stillToRead = 0;
      } 
    else
      {
      result += this->Read( buffer, sizeof(char), Self::SeekBlockSize );
      stillToRead -= Self::SeekBlockSize;
      }
    }

  this->m_BytesRead += result;
  return this->m_BytesRead;
}

