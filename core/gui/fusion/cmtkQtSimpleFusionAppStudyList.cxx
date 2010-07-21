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

#include "cmtkQtSimpleFusionApp.h"

#include <qmessagebox.h>

namespace
cmtk
{

void 
QtSimpleFusionApp::slotSetStudyList( StudyList::SmartPtr& studyList )
{
  if ( this->m_StudyList != studyList ) 
    {
    if ( studyList->size() < 2 ) 
      {
      QMessageBox::critical( NULL, "Error", "Studylist must contain at least two studies.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
      return;
      }
    
    this->m_StudyList = studyList;
    
    ReferenceStudy = this->m_StudyList->GetStudy( 0 );
    ReferenceStudy->ReadVolume( false /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
    
    MainWindow->slotDeleteAllStudies();
    StudyList::iterator it = this->m_StudyList->begin();
    while ( it != this->m_StudyList->end() )
      {
      Study::SmartPtr study = it->first;
      study->ReadVolume( false /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
      MainWindow->slotAddStudy( study );
      ++it;
      }
    }
  
  emit( signalStudyListChanged() );
}

void
QtSimpleFusionApp::slotSetReferenceStudy( Study::SmartPtr& study )
{
  ReferenceStudy = study;
  emit( signalReferenceStudyChanged() );
}

void 
QtSimpleFusionApp::slotAddStudy( Study::SmartPtr& study )
{
  if ( study ) 
    {
    study->ReadVolume( false /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
    
    if ( !study->GetVolume() ) 
      {
      QMessageBox::critical( NULL, "Error", "Cannot read image data for this study.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
      return;
      }
    
    this->m_StudyList->AddStudy( study );
    if ( this->m_StudyList->size() == 1 )
      this->ReferenceStudy = study;
    
    MainWindow->slotAddStudy( study );
    emit( signalStudyListChanged() );
    }
}

void
QtSimpleFusionApp::slotDataChanged( Study::SmartPtr& study )
{
  if ( study ) 
    {
    study->UpdateFromVolume();
    emit( signalDataChanged( study ) );
    }
}

} // namespace cmtk
