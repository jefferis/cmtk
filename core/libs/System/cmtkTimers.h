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

#ifndef __cmtkTimers_h_included_
#define __cmtkTimers_h_included_

#include <cmtkconfig.h>

#ifdef _MSC_VER
#  include <time.h>
#endif

#ifdef HAVE_TIME_H
#  include <time.h>
#endif

#ifdef HAVE_SYS_TIME_H
#  include <sys/time.h>
#endif

#ifdef HAVE_SYS_TIMES_H
#  include <sys/times.h>
#endif

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#include <iostream>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Namespace that contains functions and variables for CPU time measurements.
 */
namespace Timers 
{
  /// Get CPU time for the current process.
  extern double GetTimeProcess();

  /// Get CPU walltime for the current process.
  extern double GetWalltime();

  /** Get CPU time for the current thread.
   *\todo We need to find an equivalent implementation of this under Windows.
   */
  extern double GetTimeThread();
}  // namespace Timers

//@}

} // namespace cmtk

#endif // #ifndef __cmtkTimers_h_included_
