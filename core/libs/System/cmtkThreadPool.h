/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#ifndef __cmtkThreadPool_h_included_
#define __cmtkThreadPool_h_included_

#include <cmtkconfig.h>

#include <cmtkCannotBeCopied.h>

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
/** Class that provides a pool of continuously running threads that can be used for reducing overhead in SMP computations.
 *
 * Every instance of this class starts a pool of threads upon initialization, which can later be used to execute arbitrary tasks.
 * The threads are held via a semaphore, which are used to signal the availability of an arbitrary number of tasks. Once the
 * task semaphores have been posted, the calling function waits for the running threads to signal back the same number of
 * completed tasks via a second semaphore.
 *
 * The main advantages of this thread execution framework are:
 *
 * 1. On platforms with inefficient, slow thread creation and joining (Win32), the use and re-use of a persistent set of threads
 * greatly improves run-time efficiency and reduces overhead.
 *
 * 2. A single computation can be broken into more tasks than there are threads running, i.e., more tasks than there are CPUs
 * available. This allows for load balancing, because tasks that complete faster than others will free the executing thread, which
 * is then available to process the next available task without waiting for any other threads.
 *
 * This class provides a global thread pool, which can (and should) be shared by all computations in a process. Creating additional
 * thread pool instances should hardly ever be necessary.
 *
 * To run tasks on the global thread pool, simply create a std::vector that contains a parameter block for each task. The size of the
 * vector also determines the number of tasks to run. For example, the thread parameter could simply be a pointer to the current
 * instance ("this") of a class that acts as the client that requests a parallel computation:
 *
 *\code
 * #include <vector>
 * #include <cmtkThreadPool.h>
 *
 * class ComputationClass
 * {
 * public:
 *   typedef ComputationClass Self;
 *
 *   void ComputeUsingSMP()
 *   {
 *     // run approximately four times as many tasks as there are threads (only one task for single thread)
 *     const size_t numberOfTasks = 4 * ThreadPool::GlobalThreadPool.GetNumberOfThreads() - 3;
 *
 *     std::vector<Self*> taskParamaters( numberOfTasks, this );
 *     cmtk::ThreadPool::GlobalThreadPool.Run( ComputeTask, taskParameters );
 *   }
 *
 *   void ComputeTask( void *const arg, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t threadCnt )
 *   {
 *     Self* Caller = static_cast<Self*>( arg );
 *     // more things to do for "Caller"
 *     // taskIdx is the index of this task; taskCnt is the total number of tasks. These two determine what part of the total work must be done.
 *     // threadIdx is the index of the "physical" thread out of threadCnt threads that are running in this pool. If temporary memory is allocated
 *     // for this function, then threadIdx can be used to index this temporary storage, thus allowing us to get by with threadCnt many spaces,
 *     // rather than taskCnt many, which is usuallu much larger.
 *   }
 * };
 *\endcode
 */
class ThreadPool :
  /// Make class uncopyable via inheritance.
  private CannotBeCopied
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
	    std::vector<TParam>& taskParameters, //!< Vector of task parameter blocks, one per task.
	    const size_t numberOfTasksOverride = 0 //!< This can be used to run a smaller number of tasks than taskParameters.size(), which is useful to allow re-use of larger, allocated vector.
    );

  /// This function is run as a thread.
  void ThreadFunction( const size_t threadIdx /**!< Index of the actual thread in the pool. */ );
  
  /** Get reference to global thread pool.
   * This is shared by all functions in the process and allows re-use of the same "physical" threads 
   * for all types of computations. The thread pool itself is a local static instance within this
   * function, thus making sure it is initialized properly (see Effective C++, 3rd, Item 4).
   */
  static Self& GetGlobalThreadPool();

  /// Thread function arguments: identify pool and index of thread in it.
  class ThreadPoolArg
  {
  public:
    /// The thread pool.
    ThreadPool* m_Pool;

    /// Index of thread in pool.
    size_t m_Index;
  };

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

  /// Thread function parameters.
  std::vector<Self::ThreadPoolArg> m_ThreadArgs;
  
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
  
  /// Flag whether threads for this pool are running.
  bool m_ThreadsRunning;

  /** Flag whether threads should continue or terminate.
   * When the thread pool is destructed, this is set to "true" to wind down all running
   * threads gracefully.
   */
  bool m_ContinueThreads;

  /// Start threads for this pool.
  void StartThreads();

  /// End threads for this pool.
  void EndThreads();
};

} // namespace cmtk

/// This is the actual low-level thread function. It calls ThreadFunction() for the cmtk::ThreadPool instance given as the function parameter.
extern "C" CMTK_THREAD_RETURN_TYPE cmtkThreadPoolThreadFunction( CMTK_THREAD_ARG_TYPE arg /**!< This is a pointer to the cmtk::ThreadPool instance.*/ );

#include <cmtkThreadPool.txx>

#endif // #ifndef __cmtkThreadPool_h_included_
