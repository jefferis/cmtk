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

#include <cmtkQtProgress.h>

#include <qapplication.h>
#include <q3progressbar.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

void
QtProgress
::BeginVirtual( const float start, const float end, const float increment, const std::string& taskName )
{
  if ( this->IsTopLevel() )
    {
    if ( ProgressBar ) 
      {
      ProgressBar->setTotalSteps( 100 );
      ProgressBar->show();
      }
    
    if ( ! ProgressDialog )
      ProgressDialog = new Q3ProgressDialog( taskName.c_str(), "Cancel", 100, ParentWindow, Qt::Popup );
    
    ProgressDialog->setModal( true );
    ProgressDialog->setCaption( "Please wait" );
    ProgressDialog->setMinimumDuration( 100 );
    ProgressDialog->show();
    ProgressDialog->setTotalSteps( 100 );
    
    qApp->processEvents();
    }
}

Progress::ResultEnum
QtProgress::UpdateProgress()
{
  const int percent = static_cast<int>( 100 * this->GetFractionComplete() );
  if ( ProgressBar )
    ProgressBar->setProgress( percent );
  if ( ProgressDialog )
    ProgressDialog->setProgress( percent );

  qApp->processEvents();

  Progress::ResultEnum result = Progress::OK;
  if ( ProgressDialog )
    if ( ProgressDialog->wasCanceled() )
      result = Progress::INTERRUPT;

  return result;
}

void
QtProgress::DoneVirtual()
{
  if ( this->IsTopLevel() )
    {
    if ( ProgressBar )
      ProgressBar->reset();
    
    if ( ProgressDialog )
      ProgressDialog->hide();
    }
}

} // namespace cmtk
