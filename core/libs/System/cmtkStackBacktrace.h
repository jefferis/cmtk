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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#ifndef __cmtkStackBacktrace_h_included_
#define __cmtkStackBacktrace_h_included_

#include <cmtkconfig.h>

#ifdef HAVE_SIGNAL_H
#  include <signal.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{
/** Class for printing stack backtrace in the event of a crash.
 * A static instance of this class is used to catch SEGFAULT and other
 * termination signals and prints a stack trace.
 *\author http://www.linuxjournal.com/article/6391
 */
class
StackBacktrace
{
public:
  /// This class.
  typedef StackBacktrace Self;

  /// Constructor: register signal handler.
  StackBacktrace();

  /// Set exit code used after catching SEGFAULT or other signals.
  static void SetExitCode( const int code = 1 )
  {
    Self::ExitCode = code;
  }

private:
#ifdef HAVE_SIGNAL_H
  /// Signal handler.
#ifndef _MSC_VER
  static void SignalHandler( int sig, siginfo_t *info, void *secret );
#else
  static void SignalHandler( int sig );
#endif
#endif

  /** Exit code.
   * This defaults to "1" but can be set to "0" for CTest testing.
   */
  static int ExitCode;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkStackBacktrace_h_included_
