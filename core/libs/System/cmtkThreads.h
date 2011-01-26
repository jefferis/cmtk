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

#ifndef __cmtkThreads_h_included_
#define __cmtkThreads_h_included_

#include <cmtkconfig.h>

#include <stdlib.h>
#include <stdio.h>
#include <vector>

#include <System/cmtkThreadSystemTypes.h>
#include <System/cmtkThreadParameters.h>

namespace
cmtk
{

/** \addtogroup System */
//@{
/** Thread-related utility functions and global configuration variables.
 */
namespace 
Threads 
{

/// Check environment variables that control thread creation.
void CheckEnvironment();

/// Check whether this system supports threads.
bool Available();

/// Number of threads to run in parallel.
extern int GetNumberOfThreads();

/** Set the number of threads to run in parallel.
 * If the given parameter is less than or equal to zero, the number of 
 * threads is set to equal the number of currently available CPUs. This also
 * applies when this function is called without any parameter.
 *\param numberOfThreads Number of parallel threads to run. This is a
 * maximum value. Particular tasks may execute with less parallel threads
 * than this value.
 *\param force If this flag is off (default), then the number of parallel
 * threads is limited by the number of available CPUs on the current host.
 * If the flag is true, the limitation is not entirely lifted, but the
 * number of threads is limited by the maximum number of threads that
 * can be created by a process on the current operating system.
 *\return The actual number of threads that the library was set to.
 */
extern int SetNumberOfThreads( const int numberOfThreads = 0, const bool force = false );

/// Return number of threads allowed per process on this system.
extern int GetMaxThreads();

/// Return number of processors currently online.
extern int GetNumberOfProcessors();

/// (Maximum) number of threads to run in parallel.
extern int NumberOfThreads;

/** Specialized but more hands-on thread scheduling function.
 *\param threadCall Thread function to be called.
 *\param numberOfThreads Number of parallel threads. This parameter must not
 * exceed the number of threads set using and returned by SetNumberOfThreads.
 *\param parameters Pointer to an array of parameter blocks for all parallel
 * threads.
 *\param parameterSize Size in bytes of each thread parameter block in the
 * array pointed to by the previous parameter.
 */
void RunThreads( ThreadFunction threadCall, const unsigned numberOfThreads, void *const parameters, const size_t parameterSize );

/** Generic thread scheduling function.
 * This function created a given number of parallel threads with
 * user-provided thread-specific parameters. It then waits for all threads
 * to complete before returning. This function is a convenience wrapper for
 * the four-parameter function of the same name.
 *\param threadCall Thread function to be called.
 *\param numberOfThreads Number of parallel threads. This parameter must not
 * exceed the number of threads set using and returned by SetNumberOfThreads.
 *\param parameters Pointer to an array of parameter blocks for all parallel
 * threads. This is a template type parameter, so arbitrary parameter blocks
 * can be used.
 */
template<class T> void RunThreads( ThreadFunction threadCall, const unsigned numberOfThreads, T *const parameters )
{
  Threads::RunThreads( threadCall, numberOfThreads, parameters, sizeof( T ) );
}
}

//@}

} // namespace cmtk

#endif // #ifndef __cmtkThreads_h_included_
