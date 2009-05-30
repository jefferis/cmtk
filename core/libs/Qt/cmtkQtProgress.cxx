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
QtProgress::SetTotalStepsVirtual( const unsigned int totalSteps )
{
  if ( ProgressBar ) 
    {
    ProgressBar->setTotalSteps( totalSteps );
    ProgressBar->show();
    }
  
  if ( ! ProgressDialog )
    ProgressDialog = new Q3ProgressDialog( "This may take a bit...", "Cancel", totalSteps, ParentWindow, Qt::Popup );

  ProgressDialog->setModal( true );
  ProgressDialog->setCaption( "Please wait" );
  ProgressDialog->setMinimumDuration( 1000 );
  ProgressDialog->show();
  ProgressDialog->setTotalSteps( totalSteps );
  
  qApp->processEvents();
}

ProgressResult
QtProgress::SetProgressVirtual( const unsigned int progress )
{
  if ( ProgressBar )
    ProgressBar->setProgress( progress );
  if ( ProgressDialog )
    ProgressDialog->setProgress( progress );

  qApp->processEvents();

  ProgressResult result = PROGRESS_OK;
  if ( ProgressDialog )
    if ( ProgressDialog->wasCanceled() )
      result = PROGRESS_INTERRUPT;

  return result;
}

void
QtProgress::DoneVirtual()
{
  if ( ProgressBar )
    ProgressBar->reset();

  if ( ProgressDialog )
    ProgressDialog->hide();
}

} // namespace cmtk
