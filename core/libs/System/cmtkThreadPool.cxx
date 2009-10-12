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

namespace
cmtk
{

ThreadPool::ThreadPool( const size_t nThreads )
{
  if ( ! nThreads )
    this->m_NumberOfThreads = cmtk::Threads::GetNumberOfThreads();
  else
    this->m_NumberOfThreads = nThreads;
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
