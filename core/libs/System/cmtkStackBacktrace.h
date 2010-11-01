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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkStackBacktrace_h_included_
#define __cmtkStackBacktrace_h_included_

#include <cmtkconfig.h>

#include <signal.h>

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

  /** Exit code.
   * This defaults to "1" but can be set to "0" for CTest testing.
   */
  static int ExitCode;

  /// Print current stack backtrace.
  static void PrintBacktrace( const int levels = 0 /*!< Maximum number of levels to display (default: 0 = no limit). */ );
};

//@}

} // namespace cmtk

#ifndef _MSC_VER
/// Signal handler.
extern "C" void cmtkStackBacktraceSignalHandler( int sig, siginfo_t *info, void *secret );
#else
extern "C" void cmtkStackBacktraceSignalHandler( int sig );
#endif

#endif // #ifndef __cmtkStackBacktrace_h_included_
