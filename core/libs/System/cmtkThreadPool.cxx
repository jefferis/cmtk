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

#include <cmtkConsole.h>

namespace
cmtk
{

ThreadPool::ThreadPool( const size_t nThreads )
  : m_NumberOfTasks( 0 ),
    m_NextTaskIndex( 0 ),
    m_TaskFunction( NULL )
{
  std::cerr << "ThreadPool constructor" << std::endl;
  if ( ! nThreads )
    this->m_NumberOfThreads = cmtk::Threads::GetNumberOfThreads();
  else
    this->m_NumberOfThreads = nThreads;

  this->m_ThreadID.resize( this->m_NumberOfThreads );
#ifdef _MSC_VER
  this->m_ThreadHandles.resize( this->m_NumberOfThreads );
#endif

#ifdef CMTK_BUILD_SMP  
#ifdef CMTK_USE_THREADS
  pthread_attr_t attr;
  pthread_attr_init(&attr);
  pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

  for ( size_t idx = 0; idx < this->m_NumberOfThreads; ++idx ) 
    {
    // nothing happened yet, so set status to OK
    const int status = pthread_create( &this->m_ThreadID[idx], &attr, cmtkThreadPoolThreadFunction, static_cast<CMTK_THREAD_ARG_TYPE>( this ) );
    
    if ( status ) 
      {
      StdErr.printf( "Creation of pooled thread #%d failed with status %d.\n", idx, status );
      exit( 1 );
      }
    }
  
  pthread_attr_destroy(&attr);
#elif defined(_MSC_VER)
  for ( size_t idx = 0; idx < this->m_NumberOfThreads; ++idx ) 
    {
    // nothing happened yet, so set status to OK
    int status = 0;
    
    this->m_ThreadHandles[idx] = CreateThread( NULL /*default security attributes*/, 0/*use default stack size*/, 
					       (LPTHREAD_START_ROUTINE) cmtkThreadPoolThreadFunction, 
					       static_cast<CMTK_THREAD_ARG_TYPE>( this ),  0/*use default creation flags*/, &this->m_ThreadID[threadIdx] );
    if ( this->m_ThreadHandles[idx] == NULL ) 
      {
      status = -1;
      }
    
    if ( status ) 
      {
      StdErr.printf( "Creation of pooled thread #%d failed with status %d.\n", idx, status );
      exit( 1 );
      }
    }
#endif // #ifdef CMTK_USE_THREADS
#endif // #ifdef CMTK_BUILD_SMP
}

ThreadPool::~ThreadPool()
{
#ifdef CMTK_BUILD_SMP
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
      pthread_join( this->m_ThreadID[idx], NULL );
      }
#endif
    }
#endif
#endif // #ifdef CMTK_BUILD_SMP
  std::cerr << "ThreadPool destructor" << std::endl;
}

void
ThreadPool::ThreadFunction()
{
#ifdef CMTK_BUILD_SMP
  const size_t threadIdx = this->GetMyThreadIndex();
  while ( true )
    {
    // wait for task waiting
    this->m_TaskWaitingSemaphore.Wait();

    // lock, get, increment next task index
    this->m_NextTaskIndexLock.Lock();
    const size_t taskIdx = this->m_NextTaskIndex;
    ++this->m_NextTaskIndex;
    this->m_NextTaskIndexLock.Unlock();
    
    // call task function
    this->m_TaskFunction( this->m_TaskParameters[taskIdx], taskIdx, this->m_NumberOfTasks, threadIdx, this->m_NumberOfThreads ); 
    
    // post "task done, thread waiting"
    this->m_ThreadWaitingSemaphore.Post(); 
    }
#endif // #ifdef CMTK_BUILD_SMP
}

size_t 
ThreadPool::GetMyThreadIndex() const
{
#ifdef CMTK_BUILD_SMP
#ifdef _MSC_VER
  const DWORD self = GetCurrentThreadId();
#else
  const pthread_t self = pthread_self();
#endif // #ifdef _MSC_VER
  for ( size_t idx = 0; idx < this->m_ThreadID.size(); ++idx )
    if ( self == this->m_ThreadID[idx] )
      return idx;
#endif // #ifdef CMTK_BUILD_SMP
  return -1;
}

}

CMTK_THREAD_RETURN_TYPE
cmtkThreadPoolThreadFunction( CMTK_THREAD_ARG_TYPE arg )
{
  static_cast<cmtk::ThreadPool*>( arg )->ThreadFunction();
  return CMTK_THREAD_RETURN_VALUE;
}
