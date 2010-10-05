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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkException_h_included_
#define __cmtkException_h_included_

#include <cmtkconfig.h>

#include <exception>
#include <stdlib.h>
#include <string>

#include <System/cmtkConsole.h>

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
  Exception( const std::string& errorMsg = "", const void *const fromObject = NULL ) 
  {
    this->m_ErrorMsg = errorMsg;
    this->m_FromObject = fromObject;
  }
  
  /// Virtual destructor.
  virtual ~Exception() throw() {};

  /** Return pointer to the error message.
   */
  const std::string& GetErrorMsg() const { return this->m_ErrorMsg; }

  /** Return pointer to the object that encountered the error condition.
   */
  const void* GetFromObject() const { return this->m_FromObject; }

  /// Overwrite inherited "what" member.
  virtual const char* what() const throw()
  {
    return this->m_ErrorMsg.c_str();
  }  

private:
  /// Pointer to a string describing the condition causing the exception.
  std::string m_ErrorMsg;

  /// Pointer to the object that encountered the error condition.
  const void* m_FromObject;
};

//@}

/// Output of command line exception.
Console& operator<<( Console& console, Exception e );

} // namespace cmtk

#endif // #ifndef __cmtkException_h_included_
