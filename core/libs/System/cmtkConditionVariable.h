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

#ifndef __cmtkConditionVariable_h_included_
#define __cmtkConditionVariable_h_included_

#include <cmtkconfig.h>

#include <System/cmtkMutexLock.h>

#if defined(CMTK_USE_PTHREADS)
#  include <pthread.h>
#endif // #if defined(CMTK_USE_PTHREADS)

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Condition variable for thread synchronization.
class ConditionVariable : 
  /// Inherit mutex lock implementation.
  public MutexLock
{
public:
  /// Constructor.
  ConditionVariable()
  {
#if defined(CMTK_USE_PTHREADS)
    pthread_cond_init( &this->m_ConditionVariable, NULL );
#endif
  }
  
  /// Destructor.
  ~ConditionVariable()
  {
#if defined(CMTK_USE_PTHREADS)
    pthread_cond_destroy( &this->m_ConditionVariable );
#endif
  }
  
  /// Wait for condition variable.
  void Wait()
  {
#if defined(CMTK_USE_PTHREADS)
    pthread_cond_wait( &this->m_ConditionVariable, &this->m_MutexLock );
#endif
  }

  /** Signal condition variable.
   * This will unblock at least one waiting thread.
   */
  void Signal()
  {
#if defined(CMTK_USE_PTHREADS)
    pthread_cond_signal( &this->m_ConditionVariable );
#endif
  }

  /** Broadcast signal condition variable.
   * This will unblock all waiting threads.
   */
  void Broadcast()
  {
#if defined(CMTK_USE_PTHREADS)
    pthread_cond_broadcast( &this->m_ConditionVariable );
#endif
  }

#if defined(CMTK_USE_PTHREADS)
private:
  /// The condition variable.
  pthread_cond_t m_ConditionVariable;
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkConditionVariable_h_included_
