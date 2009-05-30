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

#ifndef __cmtkThreads_h_included_
#define __cmtkThreads_h_included_

#include <cmtkconfig.h>

#include <cmtkMemory.h>

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#ifdef _OPENMP
#  include <omp.h>
#endif 

#ifndef _MSC_VER
#  if defined(CMTK_USE_THREADS)
#    include <pthread.h>
#  endif
/// Return type of a thread function.
#  define CMTK_THREAD_RETURN_TYPE void*
/// Return type of a thread function.
#  define CMTK_THREAD_ARG_TYPE void*
/// Return value of a thread function.
#  define CMTK_THREAD_RETURN_VALUE NULL
#else // _MSC_VER
#  include <windows.h>
#  define CMTK_THREAD_RETURN_TYPE DWORD
/// Return type of a thread function.
#  define CMTK_THREAD_ARG_TYPE LPVOID
#  define CMTK_THREAD_RETURN_VALUE NULL
#endif // _MSC_VER

namespace
cmtk
{

/** \addtogroup System */
//@{

#ifndef _MSC_VER
#  if defined(CMTK_USE_THREADS)
typedef pthread_t ThreadIDType;
#  else
/// Dummy definition for non-threading builds.
typedef int ThreadIDType;
#  endif
#else // _MSC_VER
typedef DWORD ThreadIDType;
#  endif // _MSC_VER

/// Type of thread function
typedef CMTK_THREAD_RETURN_TYPE (*ThreadFunction)(CMTK_THREAD_ARG_TYPE);

#ifndef CMTK_MAX_THREADS
/** Maximum number of threads supported.
 * This value determines the size of the statically allocated array of thread
 * IDs. If you need more than 256 parallel threads (default), you may use
 * -DCMTK_MAX_THREADS=<value> as a preprocessor switch when compiling this
 * library.
 */
#define CMTK_MAX_THREADS 256
#endif

/** Base class for thread parameter blocks.
 * This doesn't hurt to have, even if we're building without thread support.
 */
template<class T> 
class ThreadParameters
{
public:
  /// Template pointer to the object starting this thread.
  T* thisObject;
  /// Unique index of this thread instance among all threads.
  unsigned int ThisThreadIndex;
  /// Total number of threads created.
  unsigned int NumberOfThreads;

  /// Thread ID of a started thread.
  ThreadIDType m_ThreadID;

#ifdef _MSC_VER
  void* m_Handle;
#endif

  /// Default constructor.
  ThreadParameters() 
  {
    this->m_ThreadID = 0;
#ifdef _MSC_VER
    m_Handle = NULL;
#endif
  }
};

/** Thread-related utility functions and global configuration variables.
 */
namespace Threads 
{

/// Initializer class.
class Initializer
{
public:
  /// Constructor: check for CMTK_NUM_THREADS environment variable
  Initializer();
};

/// Static initializer instance.
extern Initializer InitializerInstance;

/// Check whether this system supports threads.
bool Available();

/// Number of threads to run in parallel.
extern int GetNumberOfThreads();

/** Set the number of threads to run in parallel.
 * If the given parameter is less than or equal to zero, the number of 
 * threads is set to equal the number of currently available CPUs. This also
 * applies when this function is called without any parameter.
 *\param numberOfThreads Number of parallel threads to run. This is a
 * maximum value. Particular tasks may execute with less parallel threads
 * than this value.
 *\param force If this flag is off (default), then the number of parallel
 * threads is limited by the number of available CPUs on the current host.
 * If the flag is true, the limitation is not entirely lifted, but the
 * number of threads is limited by the maximum number of threads that
 * can be created by a process on the current operating system.
 *@return The actual number of threads that the library was set to.
 */
extern int SetNumberOfThreads( const int numberOfThreads = 0, const bool force = false );

/// Return number of threads allowed per process on this system.
extern int GetMaxThreads();

/// Return number of processors currently online.
extern int GetNumberOfProcessors();

/// (Maximum) number of threads to run in parallel.
extern int NumberOfThreads;

/** Specialized but more hands-on thread scheduling function.
 *@param threadCall Thread function to be called.
 *@param numberOfThreads Number of parallel threads. This parameter must not
 * exceed the number of threads set using and returned by SetNumberOfThreads.
 *@param parameters Pointer to an array of parameter blocks for all parallel
 * threads.
 *@param parameterSize Size in bytes of each thread parameter block in the
 * array pointed to by the previous parameter.
 */
void RunThreads( ThreadFunction threadCall, const unsigned numberOfThreads, void *const parameters, const size_t parameterSize );

/** Generic thread scheduling function.
 * This function created a given number of parallel threads with
 * user-provided thread-specific parameters. It then waits for all threads
 * to complete before returning. This function is a convenience wrapper for
 * the four-parameter function of the same name.
 *@param threadCall Thread function to be called.
 *@param numberOfThreads Number of parallel threads. This parameter must not
 * exceed the number of threads set using and returned by SetNumberOfThreads.
 *@param parameters Pointer to an array of parameter blocks for all parallel
 * threads. This is a template type parameter, so arbitrary parameter blocks
 * can be used.
 */
template<class T> void RunThreads( ThreadFunction threadCall, const unsigned numberOfThreads, T *const parameters )
{
  Threads::RunThreads( threadCall, numberOfThreads, parameters, sizeof( T ) );
}
}

/** Array of thread parameters.
 * This array initializes the non-template type specific fields of the thread
 * parameter structure.
 */
template<class TClass,class TParam = ThreadParameters<TClass> >
class ThreadParameterArray
{
public:
  /** Constructor.
   * Allocate array and initialize generic fields.
   */
  ThreadParameterArray
  ( TClass *const thisObject, const size_t numberOfThreads )
  {
    this->m_AsynchronousThreadsRunning = false;
    this->m_NumberOfThreads = numberOfThreads;
    this->m_Ptr = Memory::AllocateArray<TParam>( numberOfThreads );
    for ( size_t i = 0; i < numberOfThreads; ++i )
      {
      this->m_Ptr[i].thisObject = thisObject;
      this->m_Ptr[i].ThisThreadIndex = i;
      this->m_Ptr[i].NumberOfThreads = numberOfThreads;
      this->m_Ptr[i].m_ThreadID = 0;
      }
  }

  /// Destructor.
  ~ThreadParameterArray()
  {
    if ( this->m_AsynchronousThreadsRunning )
      this->CancelAsynchronousThreads();
    delete[] this->m_Ptr;
  }

  /// Constant access operator.
  const TParam& operator[]( const size_t i ) const { return this->m_Ptr[i]; }

  /// Access operator.
  TParam& operator[]( const size_t i ) { return this->m_Ptr[i]; }

  /// Return pointer to array.
  TParam* GetPtr() { return this->m_Ptr; }

  /// Return constant pointer to array.
  const TParam* GetPtr() const { return this->m_Ptr; }

  /// Return number of threads.
  size_t GetNumberOfThreads() const { return this->m_NumberOfThreads; }

  /// Run thread function in parallel.
  void RunInParallel( ThreadFunction threadCall )
  {
    Threads::RunThreads( threadCall, this->GetNumberOfThreads(), this->GetPtr(), sizeof( TParam ) );
  }

  /// Run thread function in parallel without joining.
  void RunInParallelAsynchronous( ThreadFunction threadCall );

  /// Collect (join) threads previously started by RunInParallelAsynchronous.
  void JoinAsynchronousThreads();

  /// Cancel (terminate) threads previously started by RunInParallelAsynchronous.
  void CancelAsynchronousThreads();

  /// Check if a given thread is running.
  bool IsRunning( const size_t idx )
  {
    return this->m_Ptr[idx].m_ThreadID;
  }

  /// Run thread functions using a static FIFO scheduler.
  void RunInParallelFIFO(ThreadFunction threadCall, const size_t numberOfThreadsTotal, const size_t firstThreadIdx = 0 );
    
private:
  /// Store number of threads and entries in parameter array.
  size_t m_NumberOfThreads;
  
  /// Pointer to parameter block array.
  TParam* m_Ptr;
 
  /// Flag for running asynchronous threads.
  bool m_AsynchronousThreadsRunning;
};

//@}

} // namespace cmtk

#include <cmtkThreads.txx>

#endif // #ifndef __cmtkThreads_h_included_
