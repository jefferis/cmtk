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

#include "cmtkThreadSemaphore.h"

#ifdef CMTK_USE_THREADS

#ifdef _MSC_VER
#  include "cmtkThreadSemaphoreWindows.txx"
#elif defined(__APPLE__) || defined(__CYGWIN__)
#  include "cmtkThreadSemaphoreAppleIsRetarded.txx"
#else
#  include "cmtkThreadSemaphorePOSIX.txx"
#endif

#else // #ifdef CMTK_USE_THREADS

namespace
cmtk
{

ThreadSemaphore::ThreadSemaphore( const unsigned int )
{}

ThreadSemaphore::~ThreadSemaphore()
{}

void
ThreadSemaphore::Post( const unsigned int )
{}

void
ThreadSemaphore::Wait() 
{}

}

#endif // #ifdef CMTK_USE_THREADS
