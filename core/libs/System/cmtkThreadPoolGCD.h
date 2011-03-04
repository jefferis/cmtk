/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#ifndef __cmtkThreadPoolGCD_h_included_
#define __cmtkThreadPoolGCD_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>
#include <System/cmtkSmartPtr.h>

#include <vector>

#include <dispatch/dispatch.h>

namespace
cmtk
{

/** \addtogroup System */
//@{
/** Class that provides access to Grand Central Dispatch through a thread pool-like API.
 */
class ThreadPoolGCD :
  /// Make class uncopyable via inheritance.
  private CannotBeCopied
{
public:
  /// This class.
  typedef ThreadPoolGCD Self;

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
  typedef void (*TaskFunction)( void *const args /*!< Pointer to parameter block for this task.*/,
				const size_t taskIdx /*!< Index of this task.*/,
				const size_t taskCnt /*!< Number of tasks.*/,
				const size_t threadIdx /*!< Index of the thread that is running this task.*/,
				const size_t threadCnt /*!< Number of threads in this pool.*/ );
  
  /** Constructor: create a pool of nThreads FIFO queues.
   */
  ThreadPoolGCD( const size_t nThreads = 0 /*!< Number of FIFO queues. This is the maximum number of simultaneously running tasks. */ );
  
  /// Destructor.
  ~ThreadPoolGCD();

  /// Return number of threads in the pool.
  size_t GetNumberOfThreads() const
  {
    return this->m_NumberOfThreads;
  }

  /// Run actual worker functions through running threads.
  template<class TParam> 
  void Run( Self::TaskFunction taskFunction /*!< Pointer to task function.*/,
	    std::vector<TParam>& taskParameters /*!< Vector of task parameter blocks, one per task.*/,
	    const size_t numberOfTasksOverride = 0 /*!< This can be used to run a smaller number of tasks than taskParameters.size(), which is useful to allow re-use of larger, allocated vector.*/ );
  
  /** Get reference to global thread pool.
   * This is shared by all functions in the process and allows re-use of the same "physical" threads 
   * for all types of computations. The thread pool itself is a local static instance within this
   * function, thus making sure it is initialized properly (see Effective C++, 3rd, Item 4).
   */
  static Self& GetGlobalThreadPoolGCD();

private:
  /// Number of running threads.
  size_t m_NumberOfThreads;

  /// Thread handles.
  std::vector<dispatch_queue_t> m_Queues;
};

} // namespace cmtk

#include "cmtkThreadPoolGCD.txx"

#endif // #ifndef __cmtkThreadPoolGCD_h_included_
