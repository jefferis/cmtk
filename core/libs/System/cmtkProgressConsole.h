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
 * This class displays the programm progress on the console, using cmtk::StdErr.
 * If this process is being run from inside Slicer3, output is based on Slicer's
 * XML-type progress reporting instead, and this is written to std::cout.
 *
 *\see http://www.slicer.org/slicerWiki/index.php/Slicer3:Execution_Model_Documentation#Showing_Progress_in_an_Application
 */
class ProgressConsole :
  /// Inherit generic progress indicator interface.
  public Progress
{
public:
  /// This class.
  typedef ProgressConsole Self;

  /// Superclass.
  typedef Progress Superclass;
  
  /// Default constructor: connect to progress indicator.
  ProgressConsole( const std::string& programName = std::string("") );

  /// Destructor: finish things up.
  virtual ~ProgressConsole();

  /// Output progress to console.
  virtual ResultEnum UpdateProgress();

protected:
  /// Begin a new level of progress reporting.
  virtual void BeginVirtual( const float start, const float end, const float increment, const std::string& taskName );

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
