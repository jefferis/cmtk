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

#include <cmtkThreads.h>

// GJ: use sysctl library to query cpu number on apple macosx
#ifdef __APPLE__
#include <sys/sysctl.h>
#endif

#include <stdlib.h>

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef HAVE_LIMITS_H
#  include <limits.h>
#endif

#include <iostream>

#ifdef CMTK_USE_THREADS
#  include <pthread.h>
#  include <errno.h>
#endif

#ifdef _OPENMP
#  include <omp.h>
#endif // _OPENMP

#include <algorithm>

#ifndef _SC_NPROCESSORS_ONLN
#  ifdef _SC_NPROC_ONLN
#    define _SC_NPROCESSORS_ONLN _SC_NPROC_ONLN
#  endif
#endif

#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

int Threads::NumberOfThreads = 0;

int
Threads::GetNumberOfThreads()
{
  if ( !Threads::NumberOfThreads ) 
    Threads::CheckEnvironment();
  
  if ( Threads::NumberOfThreads ) 
    {
    return Threads::NumberOfThreads;
    } 
  else
    {
    return std::min( Threads::GetNumberOfProcessors(), Threads::GetMaxThreads() );
    }
}

bool
Threads::Available()
{
#ifndef CMTK_BUILD_SMP
  // Threads are not activated in this library.
  return false;
#else

#ifdef _MSC_VER
  // we know that threads are available on Windows32.
  return true;
#elif defined (__APPLE__)
	// GJ Apple doesn't seem to provide these _SC_THREADS etc constants 
	// in unistd.h -  I don't know how you get the info dynamically
	// so just get it statically at compile time
#ifdef _POSIX_THREADS
		//#warning POSIX THREADS AVAILABLE
  return true;
#else
  return false;
#endif
#else
  return (sysconf(_SC_THREADS) == -1);
#endif
#endif
}


int
Threads
::SetNumberOfThreads
( const int numberOfThreads, const bool force )
{
  if ( numberOfThreads )
    {
    if ( force )
      {
      NumberOfThreads = std::min( numberOfThreads, Threads::GetMaxThreads() );
      }
    else
      {
      NumberOfThreads = std::min( numberOfThreads, Threads::GetNumberOfProcessors() );
      }
    }
  else
    NumberOfThreads = std::min( Threads::GetNumberOfProcessors(), Threads::GetMaxThreads() );

#ifdef _OPENMP
  omp_set_num_threads( NumberOfThreads );
#endif

  return NumberOfThreads;
}

int
Threads::GetMaxThreads()
{
#ifdef _MSC_VER
  return CMTK_MAX_THREADS;
#  elif defined(_POSIX_THREAD_THREADS_MAX)
  return  _POSIX_THREAD_THREADS_MAX;
#  elif defined(PTHREAD_THREADS_MAX)
  return  PTHREAD_THREADS_MAX;
// GJ added extra elif for Apple
// not sure if this max is sensible, but can't find
// anywhere to query it
#elif defined(__APPLE__) && defined(CMTK_USE_THREADS)
  return std::min( CMTK_MAX_THREADS, CMTK_MAX_THREADS ); 
#  elif defined(CMTK_USE_THREADS)
  const long sysconfNumThreads = sysconf( _SC_THREAD_THREADS_MAX );
  if ( sysconfNumThreads == -1 )
    return CMTK_MAX_THREADS;
  else
    return std::min( (int)sysconfNumThreads, CMTK_MAX_THREADS );
#  else
  return 1;  
#  endif
}

int
Threads::GetNumberOfProcessors()
{
#ifdef _MSC_VER 
  SYSTEM_INFO systemInfo;
  GetSystemInfo( &systemInfo ); 
  return std::min<int>( systemInfo.dwNumberOfProcessors, CMTK_MAX_THREADS );
#elif defined(__APPLE__)
  // use sysctl to get number of available cpus on apple.  Copied from:
  // developer.apple.com/documentation/Porting/Conceptual/PortingUnix/index.html
  char *name="hw.activecpu";
  int nproc, retval;  size_t len=4;
  retval=sysctlbyname(name, &nproc, &len, NULL, 0);
  return nproc;
#elif defined(CMTK_USE_THREADS)
  return sysconf( _SC_NPROCESSORS_ONLN );
#else
  return 1;
#endif
}

void
Threads::RunThreads
( ThreadFunction threadCall, const unsigned numberOfThreads, void *const parameters, const size_t parameterSize )
{
#ifdef _OPENMP
  omp_set_num_threads( std::max<int>( 1, 1+GetNumberOfThreads()-numberOfThreads ) );
#endif

#ifdef CMTK_BUILD_SMP
  ThreadIDType Thread[CMTK_MAX_THREADS];
#ifdef _MSC_VER
  HANDLE ThreadHandles[CMTK_MAX_THREADS];
#endif
#endif

#ifdef CMTK_USE_THREADS
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif

  for ( unsigned threadIdx = 1; threadIdx < numberOfThreads; ++threadIdx ) 
    {
    void *threadParameters = ((char*) parameters) + threadIdx * parameterSize;

    // nothing happened yet, so set status to OK
    int status = 0;

    // We do NOT run the last thread as a thread, but rather continue the
    // main thread. This should allow us to compute the actual real time
    // spent on the longest thread, at least approximately.
#ifdef CMTK_BUILD_SMP
#ifdef _MSC_VER
    ThreadHandles[threadIdx] = CreateThread( NULL /*default security attributes*/, 0/*use default stack size*/, (LPTHREAD_START_ROUTINE) threadCall, threadParameters,  0/*use default creation flags*/, &Thread[threadIdx] );
    if ( ThreadHandles[threadIdx] == NULL ) 
      {
      status = -1;
      }
#else // _MSC_VER
    status = pthread_create( &Thread[threadIdx], &attr, threadCall, threadParameters );
#endif // _MSC_VER
#else
    // we're not actually using SMP, so simply run everything "by hand".
    threadCall( threadParameters );
    // this should never fail really
#endif
    
    if ( status ) 
      {
      fprintf( stderr, "Creation of thread #%d failed with status %d.\n", threadIdx, status );
#if defined(CMTK_BUILD_SMP) && defined(CMTK_USE_THREADS)
      Thread[threadIdx] = 0;
#endif
      threadCall( threadParameters );
      }
    }
  
  // Run thread #0
  threadCall( parameters );
  
  // Collect thread results in reverse order. Threads with higher numbers are
  // likely to have to do an iteration less than those with lower numbers.
  // So we collect the shorter running ones first and give more time to the
  // longer running ones while dealing with the administrational overhead of
  // the ones that are already finished.
  for ( int threadIdx = numberOfThreads-1; threadIdx; --threadIdx ) 
    {
#ifdef CMTK_BUILD_SMP
#  ifdef _MSC_VER
    WaitForSingleObject( ThreadHandles[threadIdx], INFINITE /*no timeout*/ );
#  else // _MSC_VER
    void *resultThread;
    if ( Thread[threadIdx] ) 
      {
      pthread_join( Thread[threadIdx], &resultThread );
      }
#  endif // _MSC_VER
#endif
    }
  
#ifdef CMTK_BUILD_SMP
#ifndef _MSC_VER
  pthread_attr_destroy(&attr);
#endif
#endif

#ifdef _OPENMP
  omp_set_num_threads( GetNumberOfThreads() );
#endif
}

void
Threads::CheckEnvironment()
{
  const char *env = getenv( "CMTK_NUM_THREADS" );
  // check legacy variable
  if ( ! env )
    env = getenv( "IGS_NUM_THREADS" );

  if ( env )
    {
    const int numThreads = atoi( env );
    if ( numThreads )
      {
      SetNumberOfThreads( numThreads );
      // cannot use StdErr here, because it may not be initialized yet
      std::cerr << "INFO: number of threads set to " << numThreads << " according to environment variable CMTK_NUM_THREADS\n";
      }
    else
      {
      // cannot use StdErr here, because it may not be initialized yet
      std::cerr << "WARNING: environment variable CMTK_NUM_THREADS is set but does not seem to contain a number larger than 0.\n";
      }
    }
}


//@}

} // namespace cmtk
