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
//  $Revision: 583 $
//
//  $LastChangedDate: 2009-10-13 13:22:25 -0700 (Tue, 13 Oct 2009) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

namespace
cmtk
{

ThreadSemaphore::ThreadSemaphore( const unsigned int initial )
{
  this->m_Counter = initial;
  
  pthread_mutex_init( &this->m_Mutex, NULL );
  pthread_cond_init( &this->m_Condition, NULL );
}

ThreadSemaphore::~ThreadSemaphore()
{
  pthread_mutex_destroy( &this->m_Mutex );
  pthread_cond_destroy( &this->m_Condition );
}

void
ThreadSemaphore::Post( const unsigned int increment ) 
{
  pthread_mutex_lock( &this->m_Mutex );
  this->m_Counter += increment;
  pthread_mutex_unlock( &this->m_Mutex );
  
  pthread_cond_broadcast( &this->m_Condition );
}

void
ThreadSemaphore::Wait() 
{
  pthread_mutex_lock( &this->m_Mutex );
  
  while ( !this->m_Counter ) 
    pthread_cond_wait( &this->m_Condition, &this->m_Mutex );
  
  --this->m_Counter;
  pthread_mutex_unlock( &this->m_Mutex );    
}

} // namespace cmtk
