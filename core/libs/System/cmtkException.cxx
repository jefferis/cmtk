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

#include "System/cmtkException.h"

#include <cstdarg>
#include <cstdio>

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

  va_list args;
  va_start(args, errorMsgFormat);
  vsnprintf( ExceptionErrorMessageBuffer, MESSAGE_BUFFER_SIZE, errorMsgFormat, args );
  va_end(args);

  ErrorMsg = ExceptionErrorMessageBuffer;
}

char* Exception::FormatErrorMsg( const char* format, ... )
{
  va_list args;
  va_start(args, format);
  vsnprintf( ExceptionErrorMessageBuffer, MESSAGE_BUFFER_SIZE, format, args );
  va_end(args);
  
  return ExceptionErrorMessageBuffer;
}

Console& operator<<( Console& console, Exception e )
{
  console << e.what();
  return console;
}

} // namespace cmtk
