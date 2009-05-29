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

#ifndef __cmtkException_h_included_
#define __cmtkException_h_included_

#include <cmtkconfig.h>

#include <exception>

#ifndef NULL
#define NULL 0
#endif

namespace
cmtk
{

/** \addtogroup Base */
//@{

/** Exception class.
 * Instances of this class are "thrown" in case of severe program errors.
 * They can be given an optional error message and a pointer to the object
 * that encountered the fatal condition.
 */
class Exception :
  /// Inherit from C++ standard exception.
  public std::exception
{
public:
  /** Constructor.
   *@param errorMsg An optional error message describing the condition causing
   * the exception.
   *@param fromObject An optional pointer to the object that encountered the
   * condition causing the exception.
   */
  Exception( const char *errorMsg = NULL, void *fromObject = NULL ) 
  {
    ErrorMsg = errorMsg;
    FromObject = fromObject;
  }
  
  /** Constructor.
   *@param errorMsg An optional error message describing the condition causing
   * the exception.
   *@param fromObject An optional pointer to the object that encountered the
   * condition causing the exception.
   */
  Exception( const void *fromObject, const char *errorMsgFormat, ... );

  /** Return pointer to the error message.
   */
  const char* GetErrorMsg() const { return ErrorMsg; }

  /** Return pointer to the object that encountered the error condition.
   */
  const void* GetFromObject() const { return FromObject; }

  /** Utility function: Format error message.
   * This function formats an error message using the printf() syntax and
   * a variable number of parameters. It uses an internal buffer private to
   * the Exception class, so the caller does not need to provide separate
   * storage for formatting the message.
   */
  static char* FormatErrorMsg( const char* format, ... );

  /// Overwrite inherited "what" member.
  virtual const char* what() const throw()
  {
    return ErrorMsg;
  }  

private:
  /// Pointer to a string describing the condition causing the exception.
  const char* ErrorMsg;

  /// Pointer to the object that encountered the error condition.
  const void* FromObject;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkException_h_included_
