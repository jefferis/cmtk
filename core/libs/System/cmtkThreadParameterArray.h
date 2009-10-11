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

#ifndef __cmtkThreadParameterArray_h_included_
#define __cmtkThreadParameterArray_h_included_

#include <cmtkconfig.h>

#include <cmtkThreads.h>
#include <cmtkThreadParameters.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

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

#include <cmtkThreadParameterArray.txx>

#endif // #ifndef __cmtkThreadParameterArray_h_included_
