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

#ifndef __cmtkThreadPool_h_included_
#define __cmtkThreadPool_h_included_

#include <cmtkconfig.h>

#include <cmtkThreadSystemTypes.h>
#include <cmtkThreadSemaphore.h>
#include <cmtkMutexLock.h>

namespace
cmtk
{

/** \addtogroup System */
//@{
/// Class the provides a pool of continuously running threads that can be used for reducing overhead in SMP computations.
class ThreadPool
{
public:
  /// Constructor: create a pool of nThreads running threads.
  ThreadPool( const size_t nThreads );

  /// Number of running threads.
  size_t nThreads;

  /// Semaphore to signal running threads when tasks are waiting.
  ThreadSemaphore m_TaskWaitingSemaphore;

  /// Semaphore that threads use to signal when tasks they are ready for the next task.
  ThreadSemaphore m_ThreadWaitingSemaphore;

  /// Total number of tasks to execute.
  size_t m_NumberOfTasks;

  /// Index of next available task.
  size_t m_NextTaskIndex;

  /// Lock to ensure exclusive access to the task index counter.
  MutexLock m_NextTaskIndexLock;

};

} // namespace cmtk

extern "C" CMTK_THREAD_RETURN_TYPE cmtkThreadPoolThreadFunction( CMTK_THREAD_ARG_TYPE );

#endif // #ifndef __cmtkThreadPool_h_included_
