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

#include <cmtkProgress.h>
#include <cmtkConsole.h>

#ifdef CMTK_USE_OPENMP
#  include <omp.h>
#endif

namespace
cmtk
{

/** \addtogroup System */
//@{

Progress* Progress::ProgressInstance = 0;
unsigned int Progress::TotalSteps = 0;
std::string Progress::m_CurrentTaskName( "" );

void
Progress::SetTotalSteps( const unsigned int totalSteps, const std::string& currentTaskName )
{
  Self::m_CurrentTaskName = currentTaskName;
  Self::TotalSteps = totalSteps;
  if ( ProgressInstance )
    ProgressInstance->SetTotalStepsVirtual( totalSteps );
}

Progress::ResultEnum 
Progress::SetProgress( const unsigned int progress )
{
  if ( ProgressInstance )
    return ProgressInstance->SetProgressVirtual( progress );
  else
    return Self::OK;
}

void
Progress::Done()
{
  if ( ProgressInstance )
    ProgressInstance->DoneVirtual();
}

} // namespace cmtk
