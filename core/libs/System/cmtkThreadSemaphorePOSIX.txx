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
//  $Revision: 583 $
//
//  $LastChangedDate: 2009-10-13 13:22:25 -0700 (Tue, 13 Oct 2009) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#include <iostream>
#include <errno.h>
#include <stdlib.h>

namespace
cmtk
{

ThreadSemaphore::ThreadSemaphore( const unsigned int initial )
{
  if ( sem_init( &this->m_Semaphore, 0, initial ) )
    {
    std::cerr << "ERROR: sem_init failed with errno=" << errno << "\n";
    exit( 1 );
    }
}

ThreadSemaphore::~ThreadSemaphore()
{
  if ( sem_destroy( &this->m_Semaphore ) )
    {
    std::cerr << "ERROR: sem_destroy failed with errno=" << errno << "\n";
    exit( 1 );
    }
  
}

void
ThreadSemaphore::Post( const unsigned int increment )
{
  for ( unsigned int idx = 0; idx < increment; ++idx )
    {
    if ( sem_post( &this->m_Semaphore ) )
      {
      std::cerr << "ERROR: sem_post failed with errno=" << errno << "\n";
      exit( 1 );
      }
    }
}

void
ThreadSemaphore::Wait() 
{
  if ( sem_wait( &this->m_Semaphore ) )
    {
    std::cerr << "ERROR: sem_post failed with errno=" << errno << "\n";
    exit( 1 );
    }
}

}
