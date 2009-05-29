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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkStackBacktrace.h>

#ifdef HAVE_STDIO_H
#  include <stdio.h>
#endif

#ifdef HAVE_EXECINFO_H
#  include <execinfo.h>
#endif

#ifdef HAVE_UCONTEXT_H
/* get REG_EIP from ucontext.h */
#  include <ucontext.h>
#endif

#ifdef HAVE_STDLIB_H
#  include <stdlib.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

#ifdef HAVE_SIGNAL_H

#ifndef _MSC_VER
void
StackBacktrace
::SignalHandler
( int sig, siginfo_t *info, void *secret )
{
#if defined(HAVE_UCONTEXT_H) && defined(REG_EIP)
  ucontext_t *uc = (ucontext_t *)secret;
#endif
  
  /* Do something useful with siginfo_t */
  if (sig == SIGSEGV)
    {
#ifdef REG_EIP
    printf( "Caught signal %d, faulty address is %p, from %p\n", sig, info->si_addr, reinterpret_cast<void*>( uc->uc_mcontext.gregs[REG_EIP] ) );
#else
    printf( "Caught signal %d, faulty address is %p\n", sig, info->si_addr );
#endif
    }
  else
    {
    printf( "Caught signal %d\n", sig);
    }
  
#ifdef HAVE_EXECINFO_H
  void *trace[16];
  const int trace_size = backtrace(trace, 16); 
  /* overwrite sigaction with caller's address */
#ifdef REG_EIP
  trace[1] = (void *) uc->uc_mcontext.gregs[REG_EIP];
#endif
  
  char *const *messages = backtrace_symbols( trace, trace_size );
  /* skip first stack frame (points here) */
  printf("[stack] Execution path:\n");
  for ( int i=1; i<trace_size; ++i )
    printf( "[stack] %s\n", messages[i] );
#else // #ifdef HAVE_EXECINFO_H
  fputs( "Sorry, stack bracktrace() not supported on this platform.\n", stderr );
#endif // #ifdef HAVE_EXECINFO_H
  
  exit( cmtk::StackBacktrace::ExitCode );
}
#else // #ifndef _MSC_VER
void
StackBacktrace
::SignalHandler
( int sig )
{
  exit( cmtk::StackBacktrace::ExitCode );
}
#endif // #ifndef _MSC_VER
#endif

StackBacktrace
::StackBacktrace()
{
#ifdef HAVE_SIGNAL_H
#ifndef _MSC_VER
  struct sigaction sa;

  sa.sa_sigaction = StackBacktrace::SignalHandler;
  sigemptyset (&sa.sa_mask);
  sa.sa_flags = SA_RESTART | SA_SIGINFO;

  sigaction(SIGSEGV, &sa, NULL);
  sigaction(SIGUSR1, &sa, NULL);
#else
#endif // #ifndef _MSC_VER
#endif
}

int StackBacktrace::ExitCode = 1;

} // namespace cmtk
