/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#ifndef __cmtkConsole_h_included_
#define __cmtkConsole_h_included_

#include <cmtkconfig.h>

#include <System/cmtkCannotBeCopied.h>

#include <iostream>
#include <string>
#include <cassert>

#include <System/cmtkMutexLock.h>
#include <System/cmtkLockingPtr.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Standard error output console for library.
class Console :
  /// Make class uncopyable via inheritance.
  private CannotBeCopied
{
public:
  /// Constructor.
  Console( std::ostream* stream ) 
    : m_StreamP( stream )
  { 
    this->IndentLevel = 0; 
  }

  /** Get terminal line width, if possible.
   * The line width is determined using an ioctl() call, if available. Line width can be
   * set (or overridden) by the user by setting the "CMTK_CONSOLE_LINE_WIDTH" environment
   * variable to the desired number of characters per line.
   */
  size_t GetLineWidth() const;

  /// Format text with line breaks etc.
  Console& FormatText( const std::string& text /*!< The text to format with line breaks.*/,
		       const size_t margin = 0 /*!< Left margin: this many space characters are printed at the beginning of each line.*/,
		       const size_t width = 0 /*<! Line width including margin. Default: determine based on terminal width, or default to 80 if that fails.*/,
		       const int firstLine = 0 /*!< Relative indent of first line, i.e., firstLine = -margin means no indentation of first line.*/ );

  /// Formatted output.
  void printf ( const char* format, ... );

  /// Flush output stream.
  void flush() 
  { 
    if ( this->m_StreamP )
      {
      LockingPtr<std::ostream> pStream( *this->m_StreamP, this->m_MutexLock );
      pStream->flush();
      }
  }

  /// Output stream operator.
  template<class T> 
  Console& operator << ( const T data ) 
  { 
    if ( this->m_StreamP )
      {
      LockingPtr<std::ostream> pStream( *this->m_StreamP, this->m_MutexLock );
      *pStream << data; 
      }
    return *this; 
  }

  /// Increment indentation level.
  void indent() { IndentLevel += 2; }

  /// Decrement indentation level.
  void unindent() { IndentLevel -= 2; }

  /// Implicit conversion to C++ ostream.
  operator std::ostream&()
  {
    assert( this->m_StreamP ); // will fail for StdNull!
    return *(this->m_StreamP);
  }

private:
  /// The system stream that we're attaching to.
  std::ostream* m_StreamP;
  
  /// Indentiation level.
  size_t IndentLevel;

  /// Perform indentation.
  void Indent( size_t level = 0 );
  
  /// Mutex lock for thread safety.
  MutexLock m_MutexLock;
};

/// Standard error output for the library.
extern Console StdErr;

/// Standard output for the library.
extern Console StdOut;

/// Standard output for the library.
extern Console StdNull;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkConsole_h_included_
