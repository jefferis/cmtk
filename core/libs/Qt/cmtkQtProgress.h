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

#ifndef __cmtkQtProgress_h_included_
#define __cmtkQtProgress_h_included_

#include <cmtkconfig.h>

#include <cmtkProgress.h>

#include <qwidget.h>
#include <qstatusbar.h>
#include <q3progressbar.h>
#include <q3progressdialog.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/** Class for interface of progress meter to Qt.
 */
class QtProgress :
  /// Inherit from internal progress class.
  public Progress
{
public:
  /// Constructor.
  QtProgress( QWidget *const parentWindow ) {
    ParentWindow = parentWindow;
    ProgressBar = NULL;
    ProgressDialog = NULL;
    this->m_ProgressWidgetMode = PROGRESS_DIALOG;
  }

  /// Set the embedded progress bar.
  void SetProgressBar( Q3ProgressBar *const progressBar ) 
  {
    ProgressBar = progressBar;
  }
  
  /// This member function initialises the Qt progress indicator.
  virtual void BeginVirtual( const float start, const float end, const float increment, const std::string& taskName = std::string("") );

  /// This member function sets the Qt progress indicator.
  virtual Progress::ResultEnum UpdateProgress();

  /// This member function deletes the Qt progress indicator.
  virtual void DoneVirtual();

  /// Progress indicator mode.
  typedef enum { 
    PROGRESS_DIALOG,
    PROGRESS_BAR
  } ProgressWidgetMode;

  /// Set progress indicator mode.
  void SetProgressWidgetMode( ProgressWidgetMode progressWidgetMode ) 
  {
    this->m_ProgressWidgetMode = progressWidgetMode;
  }

private:
  /// The progress window parent widget.
  QWidget* ParentWindow;

  /// The progress bar widget.
  Q3ProgressBar* ProgressBar;

  /// Progress dialog widget.
  Q3ProgressDialog* ProgressDialog;

  /// Progress widget mode.
  ProgressWidgetMode m_ProgressWidgetMode;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtProgress_h_included_
