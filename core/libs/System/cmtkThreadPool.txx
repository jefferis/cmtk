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

template<class TParam> 
void
cmtk::ThreadPool::RunWorkerFunction
( const size_t numberOfTasks, TParam* taskParameters, TaskFunction taskFunction )
{
#ifdef CMTK_BUILD_SMP
  // set task function
  this->m_TaskFunction = taskFunction;

  // initialize task index and count
  this->m_NumberOfTasks = numberOfTasks;
  this->m_NextTaskIndex = 0;

  // set parameter pointers and post semaphores for tasks waiting
  for ( size_t idx = 0; idx < numberOfTasks; ++idx )
    {
    this->m_TaskParameters[idx] = &(taskParameters[idx]);
    this->m_TaskWaitingSemaphore.Post();
    }

  // now wait for all tasks to complete, as signaled via the "thread waiting" semaphore.
  for ( size_t idx = 0; idx < numberOfTasks; ++idx )
    {
    this->m_ThreadWaitingSemaphore.Wait();
    }  
#else
  // without SMP, just run everything sequentially.
  for ( size_t idx = 0; idx < numberOfTasks; ++idx )
    {
    taskFunction( taskParameters[idx], idx, numberOfTasks );
    }
#endif
}
