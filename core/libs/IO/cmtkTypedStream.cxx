/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkTypedStream.h"

#include <System/cmtkFileUtils.h>
#include <System/cmtkConsole.h>

#include <string.h>
#include <stdlib.h>
#include <cstdarg>
#include <limits.h>

#ifdef HAVE_MALLOC_H
#  include <malloc.h>
#endif

#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif

#ifdef HAVE_SYS_STAT_H
#  include <sys/stat.h>
#endif

namespace
cmtk
{

/** \addtogroup IO */
//@{

TypedStream::TypedStream()
{
  File = NULL;
  GzFile = NULL;

  PrecisionFloat = 6;
  PrecisionDouble = 10;

  this->m_Status = Self::ERROR_NONE;
  this->m_DebugFlag = Self::DEBUG_OFF;

  memset( this->Buffer, 0, sizeof( this->Buffer ) );
  this->BufferKey = this->BufferValue = NULL;

  SplitPosition = NULL;
}

int
TypedStream
::StringCmp( const char *s1, const char *s2 )
{
  for (; *s1 && *s2; s1++, s2++) 
    {
    if (*s1 == ' ' || *s1 == '\t' || *s1 == '\n' || *s2 == ' ' || *s2 == '\t' || *s2 == '\n') 
      {
      break;
      }
    if (*s1 == *s2)
      continue;
    if (*s1 >= 'a' && *s1 <= 'z') 
      {
      if (*s1 - ('a'-'A') == *s2)
	continue;
      }
    if (*s2 >= 'a' && *s2 <= 'z') 
      {
      if (*s2 - ('a'-'A') == *s1)
	continue;
      }
    return 1;
    }
  
  if ((*s1 == ' ' || *s1 == '\0' || *s1 == '\t' || *s1 == '\n') && (*s2 == ' ' || *s2 == '\0' || *s2 == '\t' || *s2 == '\n')) 
    {
    return 0;
    }
  
  return 1;
}

char*
TypedStream
::StringSplit( char * s1 ) const
{
  if (s1)
    SplitPosition = s1-1;
  if (SplitPosition == NULL)
    return NULL;
  
  /* skip over leading white space */
  for ( SplitPosition++; *SplitPosition == '\0' || *SplitPosition == ' ' || 
	  *SplitPosition == '\t' || *SplitPosition == '\n'; SplitPosition++ )
    if ( *SplitPosition == '\0' )
      return NULL;
  
  s1 = SplitPosition;
  
  /* find token's end */
  if ( *SplitPosition == '\"' ) 
    {
    /* skip over the special string token */
    for ( SplitPosition++; *SplitPosition && *SplitPosition != '\n' && *SplitPosition != '\t'; SplitPosition++) 
      {
      if ( *SplitPosition == '\\' && *(SplitPosition+1) ) 
	{
	SplitPosition++;
	continue;
	}
      if ( *SplitPosition == '\"' ) 
	{
	SplitPosition++;
	break;
	}
      }
    } 
  else
    {
    /* skip over a numeric value */
    for ( ; *SplitPosition; SplitPosition++ ) 
      {
      if ( *SplitPosition == ' ' || *SplitPosition == '\t' || *SplitPosition == '\n')
	break;
      }
    }
  
  if ( *SplitPosition ) 
    {
    *SplitPosition = '\0';
    } 
  else
    {
    SplitPosition = NULL;
    }
  
  return s1;
}

void
TypedStream
::DebugOutput( const char* format, ... )
{
  if ( this->m_DebugFlag != Self::DEBUG_ON ) return;

  static char buffer[1024];

  va_list args;
  va_start(args, format);
  vsnprintf( buffer, sizeof( buffer ), format, args );
  va_end(args);

  fputs( buffer, stderr );
  fputs( "\n", stderr );
}

} // namespace cmtk
