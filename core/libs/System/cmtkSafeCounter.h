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

#ifndef __cmtkSafeCounter_h_included_
#define __cmtkSafeCounter_h_included_

#include <cmtkconfig.h>

#include <System/cmtkMutexLock.h>
#include <System/cmtkLockingPtr.h>
#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Thread-safe counter.
 * If the library is compiled without SMP support, then this class implements
 * only an interface to a simple unsigned int counter. With SMP support, a
 * mutex lock supplied by the MutexLock class is used to protect the
 * counter from concurrent access of multiple threads.
 */
class SafeCounter
{
public:
  /// Constructor.
  SafeCounter( const unsigned int counter = 0 ) : m_Counter( counter ) {}

  /// Retrieve counter value.
  unsigned int Get() const { return this->m_Counter; }

  /// Increment and return new counter value.
  unsigned int Increment() volatile 
  { 
    LockingPtr<unsigned int> counter( this->m_Counter, this->m_Mutex );
    return ++(*counter);
  }
  
  /// Decrement and return new counter value.
  unsigned int Decrement() volatile 
  { 
    LockingPtr<unsigned int> counter( this->m_Counter, this->m_Mutex );
    return --(*counter);
  }

private:
  /// The actual counter.
  unsigned int m_Counter;

  /// Mutex for thread-safe exclusive access to counter.
  MutexLock m_Mutex;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkSafePtr_h_included_
