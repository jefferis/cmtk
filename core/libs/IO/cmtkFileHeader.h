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

#ifndef __cmtkFileHeader_h_included_
#define __cmtkFileHeader_h_included_

#include <cmtkconfig.h>

#include <System/cmtkMemory.h>

namespace
cmtk
{

/** \addtogroup IO */
//@{

/** Access to fields in a binary file header.
 */
class FileHeader
{
public:
  /// Constructor.
  FileHeader( void *const header, const bool isBigEndian = true ) 
  {
    Header = static_cast<char*>( header );
    IsBigEndian = isBigEndian;
  }
  
  /// Retrieve field of arbitrary type from header.
  template<class T> T GetField( const size_t offset ) 
  {
    T result;
    memcpy( &result, Header+offset, sizeof( T ) );
#ifdef WORDS_BIGENDIAN
    if ( ! IsBigEndian ) result = Memory::ByteSwap( result );
#else
    if ( IsBigEndian ) result = Memory::ByteSwap( result );
#endif
    return result;
  }

  /// Store field of arbitrary type to header.
  template<class T> void StoreField( const size_t offset, T value ) 
  {
#ifdef WORDS_BIGENDIAN
    if ( (sizeof(T) > 1) && ! IsBigEndian ) value = Memory::ByteSwap( value );
#else
    if ( (sizeof(T) > 1) && IsBigEndian ) value = Memory::ByteSwap( value );
#endif
    memcpy( Header+offset, &value, sizeof( T ) );
  }
  
  /// Store null-terminated string to header.
  void StoreFieldString( const size_t offset, const char* value, const size_t maxlen = 0 ) 
  {
    if ( maxlen )
      strncpy( Header+offset, value, maxlen );
    else
      strcpy( Header+offset, value );
  }
  
  /** Compare n-character string.
   *\return The comparison result code as returned by memcmp(), i.e., the result is 0 if the given string matches
   * the respective portion of the header.
   */
  bool CompareFieldStringN( const size_t offset, const char* value, const size_t len ) const
  {
    return !memcmp( Header+offset, value, len );
  }
  
  /// Retrieve array of arbitrary type from header.
  template<class T> void GetArray
  ( T *const target, const size_t offset, const size_t length = 1 ) 
  {
    memcpy( target, Header+offset, length * sizeof( T ) );
#ifdef WORDS_BIGENDIAN
    if ( ! IsBigEndian ) 
      {
      for ( size_t i = 0; i < length; ++i ) 
	target[i] = Memory::ByteSwap( target[i] );
      }
#else
    if ( IsBigEndian ) 
      {
      for ( size_t i = 0; i < length; ++i ) 
	target[i] = Memory::ByteSwap( target[i] );
      }
#endif
  }
  
private:
  /// Pointer to the binary header.
  char *Header;

  /// Flag for endianness.
  bool IsBigEndian;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkFileHeader_h_included_
