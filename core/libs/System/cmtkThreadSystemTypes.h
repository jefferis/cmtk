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

#ifndef __cmtkThreadSystemTypes_h_included_
#define __cmtkThreadSystemTypes_h_included_

#ifndef _MSC_VER
#  if defined(CMTK_USE_THREADS)
#    include <pthread.h>
#  endif
/// Return type of a thread function.
#  define CMTK_THREAD_RETURN_TYPE void*
/// Return type of a thread function.
#  define CMTK_THREAD_ARG_TYPE void*
/// Return value of a thread function.
#  define CMTK_THREAD_RETURN_VALUE NULL
#else // _MSC_VER
#  include <windows.h>
#  define CMTK_THREAD_RETURN_TYPE DWORD
/// Return type of a thread function.
#  define CMTK_THREAD_ARG_TYPE LPVOID
#  define CMTK_THREAD_RETURN_VALUE NULL
#endif // _MSC_VER

namespace
cmtk
{

/** \addtogroup System */
//@{

#ifndef _MSC_VER
#  if defined(CMTK_USE_THREADS)
typedef pthread_t ThreadIDType;
#  else
/// Dummy definition for non-threading builds.
typedef int ThreadIDType;
#  endif
#else // _MSC_VER
typedef DWORD ThreadIDType;
#  endif // _MSC_VER

/// Type of thread function
typedef CMTK_THREAD_RETURN_TYPE (*ThreadFunction)(CMTK_THREAD_ARG_TYPE);

#ifndef CMTK_MAX_THREADS
/** Maximum number of threads supported.
 * This value determines the size of the statically allocated array of thread
 * IDs. If you need more than 256 parallel threads (default), you may use
 * -DCMTK_MAX_THREADS=<value> as a preprocessor switch when compiling this
 * library.
 */
#define CMTK_MAX_THREADS 256
#endif

//@}

} // namespace cmtk

#endif // #ifndef __cmtkThreadSystemTypes_h_included_
