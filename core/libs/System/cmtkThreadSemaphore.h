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

#ifndef __cmtkThreadSemaphore_h_included_
#define __cmtkThreadSemaphore_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>
#include <System/cmtkMutexLock.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Semaphore for thread synchronization.
 * This is implemented using two low-level mutex locks. See here:
 * http://ptspts.blogspot.com/2010/05/how-to-implement-semaphore-using.html
 * http://webhome.csc.uvic.ca/~mcheng/460/notes/gensem.pdf
 */
class ThreadSemaphore :
  /// Make class uncopyable via inheritance.
  private CannotBeCopied
{
public:
  /// Initialize semaphore.
  ThreadSemaphore( const int initial = 0 ) 
    : m_Counter( initial ) 
  {
    this->m_Delay.Lock();
  }

  /// Post semaphore.
  void Post( const int increment = 1 )
  {
    this->m_Mutex.Lock();
    this->m_Counter += increment;
    
    if ( this->m_Counter <= 0 )
      this->m_Delay.Unlock();
    else
      this->m_Mutex.Unlock();
  }
  
  /// Wait for semaphore.
  void Wait()
  {
    this->m_Mutex.Lock();
    --this->m_Counter;

    if ( this->m_Counter < 0 )
      {
      this->m_Mutex.Unlock();
      this->m_Delay.Lock();
      }
    this->m_Mutex.Unlock();
  }

private:
  /// Semaphore value.
  int m_Counter;

  /// First mutex.
  MutexLock m_Mutex;

  /// Second mutex.
  MutexLock m_Delay;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkThreadSemaphore_h_included_
