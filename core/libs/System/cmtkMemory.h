/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2004-2009 SRI International
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

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

#include <exception>
#include <iostream>
#include <typeinfo>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Memory-related helper functions.
namespace Memory
{

/// Safe allocation of C++ array with catching and output of exceptions
template<class T>
inline 
T*
AllocateArray( const size_t size )
{
  try
    {
    return new T[size];
    }
  catch ( std::exception& e )
    {
    std::cerr << "cmtk::Memory::AllocateArray caught exception '" << e.what() << "' allocating " << size << " objects of type " << typeid(T).name() << std::endl;
    throw;
    }
  return NULL;
}

/** Set (fill) memory region with given value.
 * This is the templated version of ::memset()
 */
template<class T>
inline void
Set( T* const to, const T value, const size_t length ) 
{
  for ( size_t i = 0; i < length; ++i )
    to[i] = value;
}

/** Copy memory region.
 * This is the templated version of ::memcpy()
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
T ByteSwap( const T value ) 
{
  T swap = value;
  char *cptr = reinterpret_cast<char*>( &swap );
  unsigned int j = sizeof(T)-1;
  for ( unsigned int i = 0; i < j; ++i, --j ) 
    {
    char tmp = cptr[i];
    cptr[i] = cptr[j];
    cptr[j] = tmp;
    }
  return swap;
}

/// Byte-swap arbitrary value in place.
template<class T>
void ByteSwapInPlace( T& value ) 
{
  char *cptr = reinterpret_cast<char*>( &value );
  unsigned int j = sizeof(T)-1;
  for ( unsigned int i = 0; i < j; ++i, --j ) 
    {
    char tmp = cptr[i];
    cptr[i] = cptr[j];
    cptr[j] = tmp;
    }
}

/**@name Functions for tracing memory allocation and de-allocation.
 * These functions can be used to watch the amount of memory allocated and 
 * freed by other functions. Memory holes can be identified and located.
 */
//@{

/** Get amount of memory used.
 *@return The number of bytes allocated by the process the calling function
 * is located in.
 */
inline size_t Used () 
{
#ifdef HAVE_MALLINFO
  struct mallinfo stats = mallinfo();
  return stats.uordblks + stats.usmblks;
#else
  return 0;
#endif
}

/** Print memory usage.
 *@param msg An optional message to be printed with the amount of bytes 
 * allocated by the current process. This allows location of memory holes.
 * The parameter may be omitted or given as NULL if no additional message is 
 * required.
 */
inline void Info ( const char *msg = NULL ) 
{
  int used = Used();
  if (msg )
    printf("%d bytes in use %s\n",used,msg);
  else
    printf("%d bytes in use.\n",used);
}

/** Print difference of memory usage.
 *@param before Number of bytes allocated before the inspected operation. This
 * parameter can be retrieved by al call to memused().
 *@param msg Name of the operation the memory allocation of which was 
 * inspected.
 */
inline void Diff ( const size_t before, const char *msg ) 
{
  int diff = Used()-before;
  if (diff<0)
    printf("%s freed %d bytes.\n",msg,-diff);
  else
    printf("%s allocated %d bytes.\n",msg,diff);
}

} // namespace Memory

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMemory_h_included_
