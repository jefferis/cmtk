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

#ifndef __cmtkConsole_h_included_
#define __cmtkConsole_h_included_

#include <cmtkconfig.h>

#include <iostream>
#include <string>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#endif

#include <cmtkMutexLock.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Standard error output console for library.
class Console
{
public:
  /// Constructor.
  Console() 
  { 
    this->IndentLevel = 0; 
    this->DebugLevel = 0; 
#ifdef CMTK_BUILD_MPI
    this->m_RankMPI = -1;
#endif
  }

  /// Format text with line breaks etc.
  Console& FormatText( const std::string& text, //!< The text to format with line breaks.
		       const size_t margin = 0, //!< Left margin: this many space characters are printed at the beginning of each line.
		       const size_t width = 80, //<! Line width including margin.
		       const int firstLine = 0 ); //!< Relative indent of first line, i.e., firstLine = -margin means no indentation of first line.

  /// Formatted output.
  void printf ( const char* format, ... );

  /// Flush output stream.
  void flush() 
  { 
    this->m_MutexLock.Lock();
    std::cerr.flush(); 
    this->m_MutexLock.Unlock();
  }

  /// Output stream operator.
  template<class T> 
  Console& operator << ( const T data ) 
  { 
#ifdef CMTK_BUILD_MPI
    // for now, skip the output entirely if this is not the root process.
    if ( this->m_RankMPI < 0 ) this->m_RankMPI = MPI::COMM_WORLD.Get_rank();
    if ( this->m_RankMPI ) return *this;
#endif
    this->m_MutexLock.Lock();
    std::cerr << data; 
    this->m_MutexLock.Unlock();
    return *this; 
  }

  /// Increment indentation level.
  void indent() { IndentLevel += 2; }

  /// Decrement indentation level.
  void unindent() { IndentLevel -= 2; }

private:
  /// Indentiation level.
  size_t IndentLevel;

  /// Perform indentation.
  void Indent( size_t level = 0 );
  
  /// Indentiation level.
  int DebugLevel;

  /// Mutex lock for thread safety.
  MutexLock m_MutexLock;

#ifdef CMTK_BUILD_MPI
  /// MPI process rank.
  int m_RankMPI;
#endif
};

/// Standard error output for the library.
extern Console StdErr;

//@}

} // namespace cmtk

#endif // #ifndef __cmtkConsole_h_included_
