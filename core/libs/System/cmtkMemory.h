/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#ifndef __cmtkMemory_h_included_
#define __cmtkMemory_h_included_

#include <cmtkconfig.h>

#include <stdlib.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Memory-related helper functions.
namespace Memory
{

/** Utility function: get next power of two.
   * http://en.wikipedia.org/wiki/Power_of_two#Algorithm_to_find_the_next-highest_power_of_two 
   */
size_t GetNextPowerOfTwo( size_t k );

/// Memory deallocation fucntion pointer.
typedef void (*DeallocatorFunctionPointer)( void *const );

/** Memory allocation for C-style array.
 */
class ArrayC
{ 
public:
  /** Allocator function: calls malloc()
   *\warning Objects in the array will not be instantiated.
   */
  template<class T>
  static T* Allocate( const size_t size )
  {
    return static_cast<T*>( malloc( size * sizeof( T ) ) );
  }
  
  /// Delete an array allocated using Allocate().
  template<class T>
  static void Delete( T *const array )
  {
    free( array );
  }

  /** Delete an array allocated using Allocate(), but referred to by a void*.
   * This function provides a universal signature that we can use for function pointers.
   */
  static void DeleteWrapper( void *const array )
  {
    free( array );
  }
};

/** Memory allocation for C++-style array.
 */
class ArrayCXX
{ 
public:
  /** Allocator function: calls new[]()
   *\warning Objects in the array will not be instantiated.
   */
  template<class T>
  static T* Allocate( const size_t size )
  {
    return new T[size];
  }
  
  /// Delete an array allocated using Allocate().
  template<class T>
  static void Delete( T *const array )
  {
    delete[]( array );
  }

  /** Delete an array allocated using Allocate(), but referred to by a void*.
   * This function provides a universal signature that we can use for function pointers.
   */
  template<class T>
  static void DeleteWrapper( void *const array )
  {
    delete[]( static_cast<T*>( array ) );
  }
};

/** Set (fill) memory region with given value.
 * This is the templated version of memset()
 */
template<class T>
inline void
Set( T* const to, const T value, const size_t length ) 
{
  for ( size_t i = 0; i < length; ++i )
    to[i] = value;
}

/** Copy memory region.
 * This is the templated version of memcpy()
 */
template<class T>
inline void
Copy( T* const to, const T* from, const size_t length ) 
{
  for ( size_t i = 0; i < length; ++i )
    to[i] = from[i];
}

/// Byte-swap arbitrary value.
template<class T>
T ByteSwap( T value ) 
{
  char *const cptr = reinterpret_cast<char*>( &value );
  unsigned int j = sizeof(T)-1;
  for ( unsigned int i = 0; i < j; ++i, --j ) 
    {
    char tmp = cptr[i];
    cptr[i] = cptr[j];
    cptr[j] = tmp;
    }
  return value;
}

/// Byte-swap arbitrary value in place.
template<class T>
void ByteSwapInPlace( T& value ) 
{
  char *const cptr = reinterpret_cast<char*>( &value );
  unsigned int j = sizeof(T)-1;
  for ( unsigned int i = 0; i < j; ++i, --j ) 
    {
    char tmp = cptr[i];
    cptr[i] = cptr[j];
    cptr[j] = tmp;
    }
}

/**\name Functions for tracing memory allocation and de-allocation.
 * These functions can be used to watch the amount of memory allocated and 
 * freed by other functions. Memory holes can be identified and located.
 */
//@{

/** Get amount of memory used.
 *\return The number of bytes allocated by the process the calling function
 * is located in.
 */
size_t Used();

/** Print memory usage.
 *\param msg An optional message to be printed with the amount of bytes 
 * allocated by the current process. This allows location of memory holes.
 * The parameter may be omitted or given as NULL if no additional message is 
 * required.
 */
void Info ( const char *msg = NULL );

/** Print difference of memory usage.
 *\param before Number of bytes allocated before the inspected operation. This
 * parameter can be retrieved by al call to memused().
 *\param msg Name of the operation the memory allocation of which was 
 * inspected.
 */
void Diff ( const size_t before, const char *msg );

} // namespace Memory

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMemory_h_included_
