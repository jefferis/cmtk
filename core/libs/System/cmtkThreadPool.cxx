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

#include <cmtkThreadPool.h>

#include <cmtkThreads.h>

#ifdef _OPENMP
#  include <omp.h>
#endif

#include <cmtkConsole.h>

namespace
cmtk
{

ThreadPool::ThreadPool( const size_t nThreads )
{
  if ( ! nThreads )
    this->m_NumberOfThreads = cmtk::Threads::GetNumberOfThreads();
  else
    this->m_NumberOfThreads = nThreads;

  this->m_ThreadID.resize( this->m_NumberOfThreads );
#ifdef _MSC_VER
  this->m_ThreadHandles.resize( this->m_NumberOfThreads );
#endif

#ifdef _OPENMP
  omp_set_num_threads( std::max<int>( 1, 1+Threads::GetNumberOfThreads()-this->m_NumberOfThreads ) );
#endif
  
#ifdef CMTK_USE_THREADS
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
#endif

  for ( size_t idx = 0; idx < this->m_NumberOfThreads; ++idx ) 
    {
    // nothing happened yet, so set status to OK
    int status = 0;
    
    // We do NOT run the last thread as a thread, but rather continue the
    // main thread. This should allow us to compute the actual real time
    // spent on the longest thread, at least approximately.
#ifdef CMTK_BUILD_SMP
#ifdef _MSC_VER
    this->m_ThreadHandles[idx] = CreateThread( NULL /*default security attributes*/, 0/*use default stack size*/, 
					       (LPTHREAD_START_ROUTINE) cmtkThreadPoolThreadFunction, 
					       static_cast<CMTK_THREAD_ARG_TYPE>( this ),  0/*use default creation flags*/, &this->m_ThreadID[threadIdx] );
    if ( this->m_ThreadHandles[idx] == NULL ) 
      {
      status = -1;
      }
#else // _MSC_VER
    status = pthread_create( &this->m_ThreadID[idx], &attr, cmtkThreadPoolThreadFunction, static_cast<CMTK_THREAD_ARG_TYPE>( this ) );
#endif // _MSC_VER
#else
    StdErr << "ERROR: cmtk::ThreadPool class cannot be used without SMP support\n";
    exit( 1 );
#endif
    
    if ( status ) 
      {
      StdErr.printf( "Creation of pooled thread #%d failed with status %d.\n", idx, status );
      exit( 1 );
      }
    }

#ifdef CMTK_BUILD_SMP
#ifndef _MSC_VER
  pthread_attr_destroy(&attr);
#endif
#endif
}

ThreadPool::~ThreadPool()
{
#ifdef CMTK_USE_THREADS  
  for ( size_t idx = 0; idx < this->m_NumberOfThreads; ++idx ) 
    {
#ifdef _MSC_VER
    DWORD resultThread;
    TerminateThread( this->m_ThreadHandles[idx], resultThread );
#else
    if ( this->m_ThreadID[idx] ) 
      {
      pthread_cancel( this->m_ThreadID[idx] );
      }
#endif
    }
#endif
  
#ifdef _OPENMP
  omp_set_num_threads( Threads::GetNumberOfThreads() );
#endif
}


}

CMTK_THREAD_RETURN_TYPE
cmtkThreadPoolThreadFunction( CMTK_THREAD_ARG_TYPE arg )
{
  cmtk::ThreadPool* myThreadPool = static_cast<cmtk::ThreadPool*>( arg );

  while ( true )
    {
    myThreadPool->m_TaskWaitingSemaphore.Wait();
    
    myThreadPool->m_NextTaskIndexLock.Lock();
    while ( myThreadPool->m_NextTaskIndex < myThreadPool->m_NumberOfTasks )
      {
      ++myThreadPool->m_NextTaskIndex;
      myThreadPool->m_NextTaskIndexLock.Unlock();
      
      myThreadPool->m_NextTaskIndexLock.Lock();
      }
    myThreadPool->m_NextTaskIndexLock.Unlock();
    
    myThreadPool->m_ThreadWaitingSemaphore.Post();
    }
  
  return CMTK_THREAD_RETURN_VALUE;
}
