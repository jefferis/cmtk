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

#ifndef __cmtkThreadSemaphore_h_included_
#define __cmtkThreadSemaphore_h_included_

#include <cmtkconfig.h>

#if defined(CMTK_USE_THREADS)
#    include <semaphore.h>
#elif defined(_MSC_VER)
#  include <Windows.h>
#endif

#include <iostream>
#include <stdlib.h>

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
#elif defined(_MSC_VER)
    this->m_Semaphore = CreateSemaphore( NULL /*default security attributes*/, initial, 32768 /*maximum count*/, NULL /*unnamed semaphore*/ );
    
    if ( this->m_Semaphore == NULL) 
      {
      std::cerr << "CreateSemaphore error: " << GetLastError() << std::endl;
      exit( 1 );
      }
#endif
  }

  /// Destroy semaphore object.
  ~ThreadSemaphore()
  {
#if defined(CMTK_USE_THREADS)
    sem_destroy( &this->m_Semaphore );
#elif defined(_MSC_VER)
    CloseHandle( this->m_Semaphore );
#endif
  }

  /// Post semaphore.
  void Post() 
  {
#if defined(CMTK_USE_THREADS)
    sem_post( &this->m_Semaphore );
#elif defined(_MSC_VER)
    ReleaseSemaphore( this->m_Semaphore, 1, NULL);
#endif
  }
  
  /// Wait for semaphore.
  void Wait() 
  {
#if defined(CMTK_USE_THREADS)
    sem_wait( &this->m_Semaphore );
#elif defined(_MSC_VER)
    WaitForSingleObject( this->m_Semaphore, INFINITE );
#endif
  }

  /** Try to wait for semaphore.
   * Returns true if successful (semaphore was decremented), false if not successful (semaphore unchanged).
   */
  bool TryWait() 
  {
#if defined(CMTK_USE_THREADS)
    return !sem_wait( &this->m_Semaphore );
#elif defined(_MSC_VER)
    const DWORD waitResult = WaitForSingleObject( this->m_Semaphore, 0 );
    return (waitResult != WAIT_TIMEOUT);
#else
    return false;
#endif 
  }

#if defined(CMTK_USE_THREADS)
private:
  /// Opaque system semaphore object.
  sem_t m_Semaphore;
#elif defined(_MSC_VER)
private:
  /// Opaque system semaphore object.
  HANDLE m_Semaphore;
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkThreadSemaphore_h_included_
