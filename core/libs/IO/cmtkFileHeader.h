/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkFileHeader_h_included_
#define __cmtkFileHeader_h_included_

#include <cmtkconfig.h>

#include <IO/cmtkFileConstHeader.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Read and write ccess to fields in a binary file header.
 */
class FileHeader
  : public FileConstHeader
{
public:
  /// Constructor.
  FileHeader( void *const header, const bool isBigEndian = true ) : FileConstHeader( header, isBigEndian ) 
  { 
    this->m_Header = static_cast<char*>( header );
  }

  /// Store field of arbitrary type to header.
  template<class T> void StoreField( const size_t offset, T value ) 
  {
#ifdef WORDS_BIGENDIAN
    if ( (sizeof(T) > 1) && ! this->m_IsBigEndian ) value = Memory::ByteSwap( value );
#else
    if ( (sizeof(T) > 1) && this->m_IsBigEndian ) value = Memory::ByteSwap( value );
#endif
    memcpy( this->m_Header+offset, &value, sizeof( T ) );
  }
  
  /// Store null-terminated string to header.
  void StoreFieldString( const size_t offset, const char* value, const size_t maxlen = 0 ) 
  {
    if ( maxlen )
      strncpy( this->m_Header+offset, value, maxlen );
    else
      strcpy( this->m_Header+offset, value );
  }
    
private:
  /// Pointer to the binary header.
  char *m_Header;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFileHeader_h_included_
