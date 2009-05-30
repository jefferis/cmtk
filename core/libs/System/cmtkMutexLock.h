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

#ifndef __cmtkMutexLock_h_included_
#define __cmtkMutexLock_h_included_

#include <cmtkconfig.h>

#if defined(CMTK_USE_THREADS)
#  include <pthread.h>
#endif

#ifdef _MSC_VER
#  include <Windows.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Generic mutex lock.
 * This class represents a thread-model independent wrapper for locks on data
 * that requires mutually exclusive access.
 */
class MutexLock 
{
public:
  /// Constructor: initialize low-level lock.
  MutexLock() 
  {
#if defined(CMTK_USE_THREADS)
    pthread_mutex_init( &this->m_MutexLock, NULL );
#else
#ifdef _MSC_VER
    this->m_MutexObject = CreateMutex( NULL /*default security attributes*/, FALSE /*initially not owned*/, NULL /*unnamed*/ );
    if ( ! this->m_MutexObject )
      {
      StdErr << "FATAL: cannot create mutex object\n";
      exit( 1 );
      }
#endif
#endif
  }

  /// Destructor: clean up low-level lock.
  ~MutexLock() 
  {
#if defined(CMTK_USE_THREADS)
    pthread_mutex_destroy( &this->m_MutexLock );
#else
#ifdef _MSC_VER
    CloseHandle( this->m_MutexObject );
#endif
#endif
  }

  /// Lock: if already locked, wait until unlocked, then lock.
  void Lock() 
  {
#if defined(CMTK_USE_THREADS)
    pthread_mutex_lock( &this->m_MutexLock );
#else
#ifdef _MSC_VER
    WaitForSingleObject( this->m_MutexObject, INFINITE );
#endif
#endif
  }

  /// Unlock.
  void Unlock() 
  {
#if defined(CMTK_USE_THREADS)
    pthread_mutex_unlock( &this->m_MutexLock );
#else
#ifdef _MSC_VER
    ReleaseMutex( this->m_MutexObject );.
#endif
#endif
  }

#if defined(CMTK_USE_THREADS)
protected:
   /** Low-level mutex lock for POSIX threads.
    */
   pthread_mutex_t m_MutexLock;
#else
#ifdef _MSC_VER
  HANDLE m_MutexObject;
#endif
#endif
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkMutexLock_h_included_
