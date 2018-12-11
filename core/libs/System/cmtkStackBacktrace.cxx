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

#include "cmtkStackBacktrace.h"

#include <stdlib.h>
#include <stdio.h>

#ifdef HAVE_EXECINFO_H
#  include <execinfo.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

StackBacktrace
::StackBacktrace()
{
#ifndef _MSC_VER
  struct sigaction sa;

  sa.sa_sigaction = cmtkStackBacktraceSignalHandler;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = SA_RESTART | SA_SIGINFO;

  sigaction(SIGSEGV, &sa, NULL);
  sigaction(SIGUSR1, &sa, NULL);
#else
#endif // #ifndef _MSC_VER
}

void
StackBacktrace
::PrintBacktrace( const int levels )
{
#ifdef HAVE_EXECINFO_H
  void *trace[16];
  const int trace_size = backtrace(trace, 16);
  
  char *const *messages = backtrace_symbols( trace, trace_size );
  /* skip first stack frame (points here) */
  printf("[stack] Execution path:\n");

  const int printLevels = levels ? levels+1 : trace_size;
  for ( int i=1; i<printLevels; ++i )
    printf( "[stack] %s\n", messages[i] );
#else // #ifdef HAVE_EXECINFO_H
  fputs( "Sorry, stack bracktrace() not supported on this platform.\n", stderr );
#endif // #ifdef HAVE_EXECINFO_H
}

int StackBacktrace::ExitCode = 1;

} // namespace cmtk

#ifndef _MSC_VER
void
cmtkStackBacktraceSignalHandler
( int sig, siginfo_t *info, void* )
{
  /* Do something useful with siginfo_t */
  if (sig == SIGSEGV)
    {
    printf( "Caught signal %d, faulty address is %p\n", sig, info->si_addr );
    }
  else
    {
    printf( "Caught signal %d\n", sig);
    }
  
  cmtk::StackBacktrace::PrintBacktrace();
  
  exit( cmtk::StackBacktrace::ExitCode );
}
#else // #ifndef _MSC_VER
void
cmtkStackBacktraceSignalHandler
( int sig )
{
  exit( cmtk::StackBacktrace::ExitCode );
}
#endif // #ifndef _MSC_VER

