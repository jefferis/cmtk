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

#include "cmtkConsole.h"

#ifdef HAVE_UNISTD_H
#  include <unistd.h>
#endif

#ifdef HAVE_SYS_IOCTL_H
#  include <sys/ioctl.h>
#endif

#ifdef HAVE_TERMIOS_H
#  include <termios.h>
#endif

#include <algorithm>
#include <cstdarg>
#include <stdio.h>

namespace
cmtk
{

/** \addtogroup System */
//@{
Console StdErr( &std::cerr );
Console StdOut( &std::cout );
Console StdNull( NULL );

size_t
Console::GetLineWidth() const
{
  // Allow override by user
  const char *env = getenv( "CMTK_CONSOLE_LINE_WIDTH" );
  if ( env )
    {
    const size_t width = atoi( env );
    if ( width )
      return width;
    }
  
#if defined(HAVE_SYS_IOCTL_H) && defined(TIOCGWINSZ)
  struct winsize sz;
  if ( ioctl(0, TIOCGWINSZ, &sz) >= 0 )
    return sz.ws_col;
#endif
  // default to 80
  return 80;
}

Console& 
Console::FormatText( const std::string& text, const size_t margin, const size_t width, const int firstLine )
{
  // Set current indentation to first line.
  size_t currentIndent = static_cast<size_t>( std::max<int>( 0, margin + firstLine ) );

  // set the actual width
  const size_t actualWidth = width ? width : this->GetLineWidth();

  // the effective length of a line
  size_t length = static_cast<size_t>( std::max<int>( 1, actualWidth - currentIndent ) ) - 1;

  // loop while text left to write
  std::string remaining = text;
  while ( remaining.length() > length )
    {
    // first, see if we have forced line breaks within the first "length" characters.
    size_t breakAt = remaining.find_first_of( "\n\r", 0 );

    // nothing forced in range, check for available spaces before max length
    if ( (breakAt == std::string::npos) || (breakAt >= length) )
      breakAt = remaining.find_last_of( " ", length+1 );

    // if no more spaces within available length, look for first available
    if ( breakAt == std::string::npos )
      breakAt = remaining.find_first_of( " ", length+1 );

    // if no more spaces at all, bail
    if ( breakAt == std::string::npos )
      break;

    this->Indent( currentIndent );
    (*this) << remaining.substr( 0, breakAt ) << "\n";
    remaining.erase( 0, breakAt+1 );
    currentIndent = margin;
    length = static_cast<size_t>( std::max<int>( 1, actualWidth - currentIndent ) ) - 1;
    }
  
  size_t breakAt = remaining.find_first_of( "\n\r", 0 );
  while ( breakAt != std::string::npos )
    {
    this->Indent( currentIndent );
    (*this) << remaining.substr( 0, breakAt ) << "\n";
    remaining.erase( 0, breakAt+1 );
    breakAt = remaining.find_first_of( "\n\r", 0 );
    currentIndent = margin;
    }

  this->Indent( currentIndent );
  return (*this) << remaining << "\n";
}

void 
Console::printf( const char* format, ... )
{
#ifdef CMTK_USE_MPI
    // for now, skip the output entirely if this is not the root process.
  if ( this->m_RankMPI < 0 ) this->m_RankMPI = MPI::COMM_WORLD.Get_rank();
    if ( this->m_RankMPI ) return;
#endif
  char buffer[1024];

  va_list args;
  va_start(args, format);
  vsnprintf( buffer, sizeof( buffer ), format, args );
  va_end(args);

  this->Indent();
  *this << buffer;
}

void 
Console::Indent
( size_t level )
{
  if ( ! level )
    level = this->IndentLevel;
  
  for ( size_t i = 0; i < level; ++i ) 
    (*this) << " ";
}

} // namespace cmtk
