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
//  $Revision$
//
//  $LastChangedDate$
//
//  $LastChangedBy$
//
*/

#ifndef __cmtkProgressConsole_h_included_
#define __cmtkProgressConsole_h_included_

#include <cmtkconfig.h>

#include <cmtkProgress.h>

#include <string>

/** \addtogroup System */
//@{

namespace
cmtk
{

/** Progress indicator with console output.
 */
class ProgressConsole :
  /// Inherit generic progress indicator interface.
  public Progress
{
public:
  /// Default constructor: connect to progress indicator.
  ProgressConsole( const std::string& programName = std::string("") );
  
  /// Output progress to console.
  virtual ProgressResult SetProgressVirtual( const unsigned int progress );

protected:
  /// This member function can be overriden by derived classes.
  virtual void SetTotalStepsVirtual( const unsigned int steps );

  /// Clean up console output.
  void DoneVirtual();

private:
  /// Name of this program.
  std::string m_ProgramName;

  /// Process time at start of task.
  double m_TimeAtStart;

  /// Flag that indicates whether we're running inside of Slicer3.
  bool m_InsideSlicer3;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkProgressConsole_h_included_
