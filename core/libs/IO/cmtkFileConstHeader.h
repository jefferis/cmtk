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

#ifndef __cmtkFileConstHeader_h_included_
#define __cmtkFileConstHeader_h_included_

#include <cmtkconfig.h>

#include <System/cmtkMemory.h>

#include <string.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Read-only access to fields in a binary file header.
 */
class FileConstHeader
{
public:
  /// Constructor.
  FileConstHeader( const void *const header, const bool isBigEndian = true ) 
  {
    this->m_ConstHeader = static_cast<const char*>( header );
    this->m_IsBigEndian = isBigEndian;
  }
  
  /// Retrieve field of arbitrary type from header.
  template<class T> T GetField( const size_t offset ) 
  {
    T result;
    memcpy( &result, this->m_ConstHeader+offset, sizeof( T ) );
#ifdef WORDS_BIGENDIAN
    if ( ! this->m_IsBigEndian ) result = Memory::ByteSwap( result );
#else
    if ( this->m_IsBigEndian ) result = Memory::ByteSwap( result );
#endif
    return result;
  }

  /// Get null-terminated string from header.
  char* GetFieldString( const size_t offset, char* value, const size_t maxlen = 0 ) 
  {
    if ( maxlen )
      strncpy( value, this->m_ConstHeader+offset, maxlen );
    else
      strcpy( value, this->m_ConstHeader+offset );
    return value;
  }
  
  /** Compare n-character string.
   *\return The comparison result code as returned by memcmp(), i.e., the result is 0 if the given string matches
   * the respective portion of the header.
   */
  bool CompareFieldStringN( const size_t offset, const char* value, const size_t len ) const
  {
    return !memcmp( this->m_ConstHeader+offset, value, len );
  }
  
  /// Retrieve array of arbitrary type from header.
  template<class T> void GetArray
  ( T *const target, const size_t offset, const size_t length = 1 ) 
  {
    memcpy( target, this->m_ConstHeader+offset, length * sizeof( T ) );
#ifdef WORDS_BIGENDIAN
    if ( ! this->m_IsBigEndian ) 
      {
      for ( size_t i = 0; i < length; ++i ) 
	target[i] = Memory::ByteSwap( target[i] );
      }
#else
    if ( this->m_IsBigEndian ) 
      {
      for ( size_t i = 0; i < length; ++i ) 
	target[i] = Memory::ByteSwap( target[i] );
      }
#endif
  }
  
private:
  /// Pointer to the binary header.
  const char *m_ConstHeader;

protected:
  /// Flag for endianness.
  bool m_IsBigEndian;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFileConstHeader_h_included_
