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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkThreadSemaphore_h_included_
#define __cmtkThreadSemaphore_h_included_

#include <cmtkconfig.h>

#if defined(CMTK_USE_THREADS)
#    include <semaphore.h>
#endif

#ifdef HAVE_TIME_H
#  include <time.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Semaphore for thread synchronization.
class ThreadSemaphore
{
public:
  /// Initialize semaphore.
  ThreadSemaphore( const unsigned int initial = 0 )
  {
#if defined(CMTK_USE_THREADS)
    sem_init( &this->m_Semaphore, 0, initial );
#endif
  }

  /// Destroy semaphore object.
  ~ThreadSemaphore()
  {
#if defined(CMTK_USE_THREADS)
    sem_destroy( &this->m_Semaphore );
#endif
  }

  /// Post semaphore.
  void Post() 
  {
#if defined(CMTK_USE_THREADS)
    sem_post( &this->m_Semaphore );
#endif
  }
  
  /// Wait for semaphore.
  void Wait() 
  {
#if defined(CMTK_USE_THREADS)
    sem_wait( &this->m_Semaphore );
#endif
  }

  /** Try to wait for semaphore.
   * Returns true if successful (semaphore was decremented), false if not successful (semaphore unchanged).
   */
  bool TryWait() 
  {
#if defined(CMTK_USE_THREADS)
    return !sem_wait( &this->m_Semaphore );
#else
    return false;
#endif 
  }

#if defined(CMTK_USE_THREADS)
private:
  /// Opaque system semaphore object.
  sem_t m_Semaphore;
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkThreadSemaphore_h_included_
