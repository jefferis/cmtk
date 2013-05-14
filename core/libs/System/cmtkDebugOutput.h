/*
//
//  Copyright 2011, 2013 SRI International
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

#ifndef __cmtkDebugOutput_h_included_
#define __cmtkDebugOutput_h_included_

#include <cmtkconfig.h>

#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

/// Class for debug output with different levels of detail.
class DebugOutput
{
public:
  /// This class.
  typedef DebugOutput Self;

  /// Constructor.
  DebugOutput( const int level = 0 ) : m_Level( level ) {}
  
  /// Output operator.
  template<class T>
  Console& operator<<( const T data ) const
  {
    return this->GetStream() << data;
  }

  /// Flush the appropriate stream for this output object.
  void Flush()
  {
    if ( this->m_Level > Self::GetGlobalLevel() )
      StdNull.flush();
    else
      StdOut.flush();
  }

  /// Get the appropriate stream for this output object.
  Console& GetStream() const
  {
    if ( this->m_Level > Self::GetGlobalLevel() )
      return StdNull;
    else
      return StdOut;
  }

  /// Set global debug level.
  static void SetGlobalLevel( const long int level )
  {
    Self::GetGlobalLevel() = level;
  }

  /// Increment global debug level by 1.
  static void IncGlobalLevel()
  {
    ++Self::GetGlobalLevel();
  }

  /// Get global debug level (reference to static variable).
  static int& GetGlobalLevel()
  {
    static int globalLevel = 0;
    return globalLevel;
  }

private:
  /** Level for this instance.
   * Output is suppressed if this is higher than the global level.
   */
  int m_Level;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDebugOutput_h_included_
