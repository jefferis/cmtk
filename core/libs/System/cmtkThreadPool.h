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
#include <cmtkSmartPtr.h>

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
  /// This class.
  typedef ThreadPool Self;

  /// Smart pointer.
  typedef SmartPointer<Self> SmartPtr;

  /** Task function: this is the interface for the functions called by the pooled threads to do the actual work.
   * The task function receives five parameters: a) a pointer to its parameter black, b) the index of the task,
   * c) the number of tasks, d) the index of the thread within the pool that is calling the task, and e) the
   * number of threads in the pool. Whereas the function should use b) and c) to determine what portion of work 
   * it needs to do, d) and e) must be used to determine, for example, what local memory should be used, if
   * temporary storage has been allocated for each thread. Because the number of tasks typically exceeds the
   * number of threads, this is more efficient than allocating temporary storage for each task.
   */
  typedef void (*TaskFunction)( void *const args, //!< Pointer to parameter block for this task.
				const size_t taskIdx, //!< Index of this task.
				const size_t taskCnt, //!< Number of tasks.
				const size_t threadIdx, //!< Index of the thread that is running this task.
				const size_t threadCont //!< Number of threads in this pool.
    );

  /** Constructor: create a pool of nThreads running threads.
   *\param nThreads Number of threads to create for this pool. By default, the number
   * of threads created is the current number of available threads, i.e., typically
   * the number of CPUs minus the number of currently running threads, if any.
   */
  ThreadPool( const size_t nThreads = 0 );

  /// Destructor: stop all running threads.
  ~ThreadPool();

  /// Return number of threads in the pool.
  size_t GetNumberOfThreads() const
  {
    return this->m_NumberOfThreads;
  }

  /// Run actual worker functions through running threads.
  template<class TParam> 
  void Run( TaskFunction taskFunction, //!< Pointer to task function.
	    const size_t numberOfTasks, //!< Number of tasks. Ideally, for scheduling, this is larger than the number of running threads, but this is not required.
	    TParam* taskParameters //!< Pointer to array of task parameters.
    );

  /// This function is run as a thread.
  void ThreadFunction();
  
private:
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

  /// The current task function.
  TaskFunction m_TaskFunction;

  /// Task function parameter pointers.
  std::vector<void*> m_TaskParameters;

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
  
  /// Get index of currently running thread (called from inside ThreadFunction()).
  size_t GetMyThreadIndex() const;
};

} // namespace cmtk

/// This is the actual low-level thread function. It calls ThreadFunction() for the cmtk::ThreadPool instance given as the function parameter.
extern "C" CMTK_THREAD_RETURN_TYPE cmtkThreadPoolThreadFunction( CMTK_THREAD_ARG_TYPE arg //!< This is a pointer to the cmtk::ThreadPool instance.
  );

#include <cmtkThreadPool.txx>

#endif // #ifndef __cmtkThreadPool_h_included_
