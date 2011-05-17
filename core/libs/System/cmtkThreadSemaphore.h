/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
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

#ifndef __cmtkThreadSemaphore_h_included_
#define __cmtkThreadSemaphore_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>

#if defined(CMTK_USE_THREADS)
#  if defined(__APPLE__) || defined(__CYGWIN__)
#    include <pthread.h>
#  else
#    include <semaphore.h>
#  endif
#elif defined(_MSC_VER)
#  include <Windows.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Semaphore for thread synchronization.
 * Because apparently Apple engineers are incapable of implementing an interface for unnamed
 * semaphores as provided by <semaphore.h>, we are building the semaphore ourselves on the
 * Mac OS platform using a mutex and a condition variable.
 */
class ThreadSemaphore :
  /// Make class uncopyable via inheritance.
  private CannotBeCopied
{
public:
  /// Initialize semaphore.
  ThreadSemaphore( const unsigned int initial = 0 );

  /// Destroy semaphore object.
  ~ThreadSemaphore();

  /// Post semaphore.
  void Post( const unsigned int increment = 1 );
  
  /// Wait for semaphore.
  void Wait();

#if defined(CMTK_USE_THREADS)
#  if defined(__APPLE__) || defined(__CYGWIN__)
private:
  /// Counter (Apple only).
  long int m_Counter;

  /// Counter mutex lock (Apple only).
  pthread_mutex_t m_Mutex;
  
  /// Condition variable (Apple only).
  pthread_cond_t m_Condition;
#  else // POSIX
  /// Opaque system semaphore object (POSIX only).
  sem_t m_Semaphore;
#  endif
#elif defined(_MSC_VER)
private:
  /// Opaque system semaphore object (Windows native only).
  HANDLE m_Semaphore;
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkThreadSemaphore_h_included_
