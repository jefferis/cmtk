/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//  Copyright 2009 SRI International
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

#ifndef __cmtkHashMapSTL_h_included_
#define __cmtkHashMapSTL_h_included_

#include <cmtkconfig.h>

#if defined(HAVE_UNORDERED_MAP)
#  include <unordered_map>
#elif defined(HAVE_UNORDERED_MAP_TR1)
#  include <tr1/unordered_map>
#elif defined(HAVE_HASH_MAP_H)
#  include <hash_map.h>
#elif defined(HAVE_HASH_MAP)
#  include <hash_map>
#endif

namespace
cmtk
{

/** Wrapper class for STL hash_map or unordered_map classes.
 * This class will use whatever hash map implementation is provided by the
 * currently used STL. If both unordered_map and hash_map are provided, the
 * former is used as it is the future standard class.
 */
template<class TKey, class TValue>
class HashMapSTL : 
#if defined(HAVE_UNORDERED_MAP) || defined(HAVE_UNORDERED_MAP_TR1)
    /// Inherit STL hash/unordered map.
    public std::unordered_map<TKey,TValue>
#elif defined(__GNUC__) && ! defined(__INTEL_COMPILER)
    public __gnu_cxx::hash_map
#else
    public std::hash_map
#endif
{
};
#else

#if SIZEOF_LONG != 8
#if defined(__GNUC__) && ! defined(__INTEL_COMPILER)
namespace __gnu_cxx
{

  template<>
    struct hash<unsigned long long>
    {
      size_t
      operator()(unsigned long long __x) const
      { return (__x & 0xffff) ^ (__x >> 32 ); }
    };

}
#endif // defined(__GNUC__) && ! defined(__INTEL_COMPILER)
#endif

#endif // #ifndef __cmtkHashMapSTL_h_included_
