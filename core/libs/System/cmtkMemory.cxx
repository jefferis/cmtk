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

#include "cmtkMemory.h"

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

#include <stdio.h> // thanks to Hans Johnson for pointing this out
#include <limits.h>

namespace
cmtk
{

namespace
Memory
{

size_t
GetNextPowerOfTwo( size_t k )
{

// http://en.wikipedia.org/wiki/Power_of_two#Algorithm_to_find_the_next-highest_power_of_two 

  if (k == 0)
    return 1;
  
  k--;
  for (size_t i=1; i<sizeof(size_t)*CHAR_BIT; i<<=1)
    k = k | k >> i;

  return k+1;
}

size_t
Used () 
{
#ifdef HAVE_MALLINFO
  struct mallinfo stats = mallinfo();
  return stats.uordblks + stats.usmblks;
#else
  return 0;
#endif
}

void
Info ( const char *msg ) 
{
  const int used = Used();
  if (msg )
    printf("%d bytes in use %s\n",used,msg);
  else
    printf("%d bytes in use.\n",used);
}

void
Diff ( const size_t before, const char *msg ) 
{
  const int diff = Used()-before;
  if (diff<0)
    printf("%s freed %d bytes.\n",msg,-diff);
  else
    printf("%s allocated %d bytes.\n",msg,diff);
}

} // namespace Memory

} // namespace cmtk
