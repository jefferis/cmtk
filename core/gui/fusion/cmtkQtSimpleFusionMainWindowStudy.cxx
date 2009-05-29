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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkQtSimpleFusionMainWindow.h>

#include <cmtkQtStudyWidget.h>
#include <cmtkQtSimpleFusionApp.h>

namespace
cmtk
{

void
QtSimpleFusionMainWindow::slotAddStudy( Study::SmartPtr& study )
{
  QtStudyWidget* studyWidget = new QtStudyWidget( this );
  studyWidget->slotSetStudy( study );

  QObject::connect( this->FusionApp, SIGNAL( signalDataChanged( Study::SmartPtr& ) ), studyWidget, SLOT( slotDataChanged( Study::SmartPtr& ) ) );

  const char* studyName = study->GetName();
  if ( ! studyName ) 
    studyName = study->SetMakeName();
  const char* defaultName = study->GetName();

  int collisionIdx = 1;
  bool collision = true;
  while ( collision )
    {
    collision = false;
    for ( int idx = 0; idx < StudyTabs->count(); ++idx )
      {
      QtStudyWidget* qtstudy = dynamic_cast<QtStudyWidget*>( StudyTabs->page( idx ) );
      if ( ! strcmp( studyName, qtstudy->GetStudy()->GetName() ) )
	{
	collision = true;
	break;
	}
      }
    
    if ( collision )
      {
      studyName = study->SetMakeName( defaultName, collisionIdx++ );
      }
    }

  StudyTabs->addTab( studyWidget, studyName );
  StudyTabs->showPage( studyWidget );

  CurrentStudy = study;
}

void
QtSimpleFusionMainWindow::slotDeleteStudy( Study::SmartPtr& study )
{
}

void
QtSimpleFusionMainWindow::slotDeleteAllStudies()
{
  while ( StudyTabs->count() )
    StudyTabs->removePage( StudyTabs->page( 0 ) );
}

} // namespace cmtk
