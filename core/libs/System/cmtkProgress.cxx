/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkProgress.h"
#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

Progress* Progress::ProgressInstance = 0;

void
Progress::Begin
( const double start, const double end, const double increment, const std::string& taskName ){
  if ( ProgressInstance )
    {
    ProgressInstance->BeginVirtual( start, end, increment, taskName );
    }
}

void
Progress::BeginVirtual
( const double start, const double end, const double increment, const std::string& taskName )
{
  this->m_RangeStack.push_front( Self::Range( start, end, increment, taskName ) );
}

void
Progress::SetProgressCurrent( const double progress )
{
  RangeStackType::iterator current = this->m_RangeStack.begin();
  
  if ( current != this->m_RangeStack.end() )
    {
    current->m_Current = progress;
    }
}

Progress::ResultEnum 
Progress::SetProgress( const double progress )
{
  if ( ProgressInstance )
    {
    ProgressInstance->SetProgressCurrent( progress );
    return ProgressInstance->UpdateProgress();
    }
  else
    return Self::OK;
}

void
Progress::Done()
{
  if ( ProgressInstance )
    ProgressInstance->DoneVirtual();
}

void
Progress::DoneVirtual()
{
  if ( this->m_RangeStack.begin() != this->m_RangeStack.end() )
    this->m_RangeStack.pop_front();
}

double
Progress::Range::GetFractionComplete( const double nestedFraction ) const
{
  return ( ( this->m_Current + nestedFraction * this->m_Increment ) - this->m_Start) / (this->m_End - this->m_Start);
}

double
Progress::GetFractionComplete() const
{
  double fraction = 0;
  for ( RangeStackType::const_iterator it = this->m_RangeStack.begin(); it != this->m_RangeStack.end(); ++it )
    {
    fraction = it->GetFractionComplete( fraction );
    }

  return fraction;
}

const std::string
Progress::GetCurrentTaskName() const
{
  RangeStackType::const_reverse_iterator it = this->m_RangeStack.rbegin();
  if ( it != this->m_RangeStack.rend() )
    return it->m_TaskName;
  return std::string("");
}

} // namespace cmtk
