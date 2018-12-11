/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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
//  $Revision: 2398 $
//
//  $LastChangedDate: 2010-10-05 14:54:37 -0700 (Tue, 05 Oct 2010) $
//
//  $LastChangedBy: torstenrohlfing $
//
*/

#ifndef __cmtkExitException_h_included_
#define __cmtkExitException_h_included_

#include <cmtkconfig.h>

#include <exception>
#include <stdlib.h>
#include <string>

#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/** Program exit exception class.
 * Throwing an object of this class should return control to the global main() function,
 * which should exit() with the given exit status. This should ensure that all local
 * automatic as well as global static objects have their destructors called.
 *
 * The global main() function can make use of this exception via the following wrapper
 * construct,
 *\code
 * int
 * main( const int argc, const char* argv[] )
 * {
 *   int exitCode = 0;
 *   try
 *     {
 *     exitCode = doMain( argc, argv );
 *     }
 *   catch ( const cmtk::ExitException& ex )
 *     {
 *     exitCode = ex.ExitCode();
 *     }
 *   return exitCode;
 * }
 *\endcode
 * where the actual function of the program is moved into the doMain() implementation function.
 */
class ExitException :
  /// Inherit from C++ standard exception.
  public std::exception
{
public:
  /** Constructor.
   */
  ExitException( const int exitCode = 0 /*!< Program exit code. */ ) : m_ExitCode( exitCode ) {}
  
  /// Virtual destructor.
  virtual ~ExitException() throw() {};

  /// Overwrite inherited "what" member.
  virtual const char* what() const throw()
  {
    return "cmtk::ExitException";
  }  

  /// Return exit code.
  int ExitCode() const
  {
    return this->m_ExitCode;
  }

private:
  /// Program exit code.
  const int m_ExitCode;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkExitException_h_included_
