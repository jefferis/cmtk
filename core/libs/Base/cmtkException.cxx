/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include <cmtkException.h>

#ifdef HAVE_STDARG_H
#  include <stdarg.h>
#  include <stdio.h>
#else
#  include <string.h>
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{

#define MESSAGE_BUFFER_SIZE 256

static char ExceptionErrorMessageBuffer[MESSAGE_BUFFER_SIZE];

Exception::Exception( const void *fromObject, const char *errorMsgFormat, ... )
{
  FromObject = fromObject;

#ifdef HAVE_STDARG_H
  va_list args;
  va_start(args, errorMsgFormat);
  vsnprintf( ExceptionErrorMessageBuffer, MESSAGE_BUFFER_SIZE, errorMsgFormat, args );
  va_end(args);
#else
  strncpy( ExceptionErrorMessageBuffer, errorMsgFormat, MESSAGE_BUFFER_SIZE );
  ExceptionErrorMessageBuffer[MESSAGE_BUFFER_SIZE-1] = 0;
#endif

  ErrorMsg = ExceptionErrorMessageBuffer;
}

char* Exception::FormatErrorMsg( const char* format, ... )
{
#ifdef HAVE_STDARG_H
  va_list args;
  va_start(args, format);
  vsnprintf( ExceptionErrorMessageBuffer, MESSAGE_BUFFER_SIZE, format, args );
  va_end(args);
#else
  strncpy( ExceptionErrorMessageBuffer, format, MESSAGE_BUFFER_SIZE );
  ExceptionErrorMessageBuffer[MESSAGE_BUFFER_SIZE-1] = 0;
#endif
  
  return ExceptionErrorMessageBuffer;
}

} // namespace cmtk
