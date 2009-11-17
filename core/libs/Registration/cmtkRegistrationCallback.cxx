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

#include <cmtkRegistrationCallback.h>

#ifdef HAVE_STDARG_H
#  include <stdarg.h>
#endif

#ifdef HAVE_SIGNAL_H
#  include <signal.h>
#endif

#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

/// Flag that is set upon SIGINT signal.
static bool InterruptSignalReceived;

RegistrationCallback::RegistrationCallback() 
{ 
#ifndef DEBUG
  InterruptSignalReceived = false;
#ifndef _MSC_VER
  signal( SIGINT, cmtkRegistrationCallbackDispatchSIGINT );
#endif
#endif
}

RegistrationCallback::~RegistrationCallback() 
{ 
#ifndef DEBUG
#  ifdef HAVE_SIGRELSE
  sigrelse( SIGINT );
#  endif
#endif
}

CallbackResult 
RegistrationCallback::ExecuteWithData
( const CoordinateVector&, const double )
{
  return InterruptSignalReceived ? CALLBACK_INTERRUPT : CALLBACK_OK;
}

CallbackResult
RegistrationCallback::Execute ()
{
  return InterruptSignalReceived ? CALLBACK_INTERRUPT : CALLBACK_OK;
}

void
RegistrationCallback::Comment ( const char* )
{
  return;
}

void
RegistrationCallback::FormatComment( const char* format, ... )
{
#ifdef HAVE_STDARG_H
  static char buffer[1024];

  va_list args;
  va_start(args, format);
  vsnprintf( buffer, sizeof( buffer ), format, args );
  va_end(args);

  this->Comment( buffer );
#endif
}

} // namespace cmtk

void
cmtkRegistrationCallbackDispatchSIGINT( int sig )
{
  if ( cmtk::InterruptSignalReceived )
    {
    cmtk::StdErr.printf( "Received repeated INT signal... exiting.\n" );
    exit( 3 );
    }

#ifndef DEBUG
  cmtk::InterruptSignalReceived = true;
#ifndef _MSC_VER
  signal( sig, cmtkRegistrationCallbackDispatchSIGINT );
#endif
#endif
  cmtk::StdErr.printf( "Received INT (%d) signal... preparing exit. Press Ctrl-C again to abort immediately.\n", sig );
}
