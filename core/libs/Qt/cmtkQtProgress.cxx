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

#include "Qt/cmtkQtProgress.h"

#include <qapplication.h>
#include <qprogressbar.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtProgress::QtProgress
( QWidget *const parentWindow ) 
{
  ParentWindow = parentWindow;
  ProgressBar = NULL;
  ProgressDialog = NULL;
  this->m_ProgressWidgetMode = PROGRESS_DIALOG;
}

void
QtProgress
::BeginVirtual( const double start, const double end, const double increment, const std::string& taskName )
{
  this->Superclass::BeginVirtual( start, end, increment, taskName );

  if ( this->IsTopLevel() )
    {
    if ( ProgressBar ) 
      {
      ProgressBar->setRange( 0, 100 );
      ProgressBar->show();
      }
    
    if ( ! ProgressDialog )
      ProgressDialog = new QProgressDialog( taskName.c_str(), "Cancel", 0, 100, ParentWindow, Qt::Dialog );
    
    ProgressDialog->setWindowModality(Qt::WindowModal);
    ProgressDialog->setModal( true );
    ProgressDialog->setMinimumDuration( 100 );
    ProgressDialog->show();
    ProgressDialog->setRange( 0, 100 );
    
    qApp->processEvents();
    }

  Progress::SetProgressInstance( this );
}

Progress::ResultEnum
QtProgress::UpdateProgress()
{
  const int percent = static_cast<int>( 100 * this->GetFractionComplete() );
  if ( ProgressBar )
    ProgressBar->setValue( percent );
  if ( ProgressDialog )
    ProgressDialog->setValue( percent );

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
