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

#include <vector>

namespace
cmtk
{

/** \addtogroup System */
//@{
/// Class the provides a pool of continuously running threads that can be used for reducing overhead in SMP computations.
class ThreadPool
{
public:
  /** Constructor: create a pool of nThreads running threads.
   *\param nThreads Number of threads to create for this pool. By default, the number
   * of threads created is the current number of available threads, i.e., typically
   * the number of CPUs minus the number of currently running threads, if any.
   */
  ThreadPool( const size_t nThreads = 0 );

  /// Destructor: stop all running threads.
  ~ThreadPool();

  /// Run actual worker functions through running threads.
  template<class TParam> void RunWorkerFunction( const size_t numberOfTasks, //!< Number of tasks. Ideally, for scheduling, this is larger than the number of running threads, but this is not required.
						 TParam* taskParameters, //!< Pointer to array of task parameters.
						 void (*taskFunction)( TParam *const args ) //!< Pointer to task function.
    );

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

private:
  /// Number of running threads.
  size_t m_NumberOfThreads;

#ifdef CMTK_BUILD_SMP
  /// Thread handles.
  std::vector<ThreadIDType> m_ThreadID;
  
#ifdef _MSC_VER
  /// Windows thread handles
  std::vector<HANDLE> m_ThreadHandles;
#endif
#endif  
};

} // namespace cmtk

extern "C" CMTK_THREAD_RETURN_TYPE cmtkThreadPoolThreadFunction( CMTK_THREAD_ARG_TYPE );

#include <cmtkThreadPool.txx>

#endif // #ifndef __cmtkThreadPool_h_included_
