/*
//
//  Copyright 2011 SRI International
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

  /// Set debug level.
  void SetLevel( const int level )
  {
    Self::m_Level = level;
  }

  /// Get debug level.
  int GetLevel()
  {
    return Self::m_Level;
  }

  /// Get global debug output object.
  static Self& GetGlobalStream()
  {
    static Self globalStream;
    return globalStream;
  }

private:
  /// Current toolkit-wide debug level.
  int m_Level;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkDebugOutput_h_included_
