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

#include <Windows.h>

#include <iostream>

namespace
cmtk
{

ThreadSemaphore::ThreadSemaphore( const unsigned int initial )
{
  this->m_Semaphore = CreateSemaphore( NULL /*default security attributes*/, initial, 32768 /*maximum count*/, NULL /*unnamed semaphore*/ );
  
  if ( this->m_Semaphore == NULL) 
    {
    std::cerr << "CreateSemaphore error: " << GetLastError() << std::endl;
    exit( 1 );
    }
}

ThreadSemaphore::~ThreadSemaphore()
{
  CloseHandle( this->m_Semaphore );
}

void
ThreadSemaphore::Post( const unsigned int increment ) 
{
  ReleaseSemaphore( this->m_Semaphore, increment, NULL);
}

void
ThreadSemaphore::Wait() 
{
  WaitForSingleObject( this->m_Semaphore, INFINITE );
}

} // namespace cmtk
