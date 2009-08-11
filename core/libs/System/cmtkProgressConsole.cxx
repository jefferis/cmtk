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

#include <cmtkProgressConsole.h>

#include <cmtkConsole.h>
#include <cmtkTimers.h>

#ifdef CMTK_USE_OPENMP
#  include <omp.h>
#endif

#include <cstdlib>

namespace
cmtk
{

/** \addtogroup System */
//@{

ProgressConsole::ProgressConsole( const std::string& programName )
  : m_ProgramName( programName )
{
  Progress::SetProgressInstance( this );

  this->m_InsideSlicer3 = ( getenv( "Slicer3_HOME" ) != NULL );
}

Progress::ResultEnum
ProgressConsole::SetProgressVirtual( const unsigned int progress )
{
#ifdef CMTK_USE_OPENMP
  if ( omp_get_thread_num() != 0 )
    {
    return Self::OK;
    }

  const float fraction = static_cast<float>( progress * omp_get_num_threads() ) / this->TotalSteps;
#else
  const float fraction = static_cast<float>( progress ) / this->TotalSteps;
#endif

  if ( this->m_InsideSlicer3 )
    {
    std::cout << "<filter-progress>" << fraction << "</filter-progress>\n";
    }
  else
    {
    if ( Self::m_CurrentTaskName.length() )
      {
      StdErr.printf( "%s: %d %%\r", Self::m_CurrentTaskName.c_str(), static_cast<int>( 100.0 * fraction ) );
      }
    else
      {
      StdErr.printf( "%d %%\r", static_cast<int>( 100.0 * fraction ) );
      }
    }

  return Self::OK;
}

void
ProgressConsole::SetTotalStepsVirtual( const unsigned int )
{
  this->m_TimeAtStart = Timers::GetTimeProcess();

  if ( this->m_InsideSlicer3 )
    {
    std::cout << "<filter-start>\n"
	      << "<filter-name>" << this->m_ProgramName << "</filter-name>\n"
	      << "<filter-comment> \"" << this->m_CurrentTaskName << "\" </filter-comment>\n"
	      << "</filter-start>\n";
    }
}

void
ProgressConsole::DoneVirtual()
{
  if ( this->m_InsideSlicer3 )
    {
    std::cout << "<filter-end>\n"
	      << "<filter-name>" << this->m_ProgramName << "</filter-name>\n"
	      << "<filter-time>" << Timers::GetTimeProcess() - this->m_TimeAtStart << "</filter-time>\n"
	      << "</filter-end>\n";
    }
  else
    {
    StdErr << "done.\n";
    }
}

} // namespace cmtk
