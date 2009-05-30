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

#include <cmtkQtTriplanarViewer.h>

#include <q3popupmenu.h>
#include <qmenubar.h>
#include <qapplication.h>
#include <q3filedialog.h>
#include <qlayout.h>
#include <Q3VBoxLayout>

#include <cmtkQtImageOperators.h>

#include <cmtkAnatomicalOrientation.h>
#include <cmtkStudy.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtTriplanarViewer::QtTriplanarViewer()
  : WindowLevel( NULL )
{
  this->setCaption( "Triplanar Image Viewer" );
  Q3PopupMenu* StudyMenu = new Q3PopupMenu();
  StudyMenu->insertItem( "&Load File...", this, SLOT( slotLoadFile() )  );
  StudyMenu->insertItem( "Load Stud&y...", this, SLOT( slotLoadStudy() )  );
  StudyMenu->insertItem( "&Reload Data...", this, SLOT( slotReloadData() )  );
  StudyMenu->insertSeparator();
  StudyMenu->insertItem( "&Save" );
  StudyMenu->insertItem( "Save &as..." );
  StudyMenu->insertItem( "&Export landmarks..." );
  StudyMenu->insertSeparator();
  StudyMenu->insertItem( "&Quit", qApp, SLOT( quit() ) );
  
  QtImageOperators* imageOperators = new QtImageOperators( &this->m_Study, this, NULL /*progressInstance*/ );
  QObject::connect( imageOperators, SIGNAL( dataChanged( Study::SmartPtr& ) ), this, SLOT( slotDataChanged( Study::SmartPtr& ) ) );
  
  // insert "Study" as first menu.
  MenuBar->insertItem( "&Study", StudyMenu, -1, 0 );
  // insert "Operators" after "View".
  MenuBar->insertItem( "&Operators", imageOperators->CreatePopupMenu(), -1, 2 );
  MenuBar->show();
  
  this->m_ImagesTab = new QWidget( this->m_ControlsTab );
  this->m_ControlsTab->addTab( this->m_ImagesTab, "Images" );
  this->m_ControlsTab->setTabEnabled( this->m_ImagesTab, false );

  this->m_StudiesListBox = new Q3ListBox( this->m_ImagesTab );
  this->m_StudiesListBox->setSelectionMode( Q3ListBox::Single );
  this->m_StudiesListBox->setColumnMode( Q3ListBox::FixedNumber );
  QObject::connect( this->m_StudiesListBox, SIGNAL( highlighted( const QString& ) ), this, SLOT( slotSwitchStudy( const QString& ) ) );

  Q3VBoxLayout* studiesLayout = new Q3VBoxLayout( this->m_ImagesTab );
  studiesLayout->setContentsMargins( 5, 5, 5, 5 );
  studiesLayout->setSpacing( 5 );

  studiesLayout->addWidget( this->m_StudiesListBox );

  QPushButton* copyColormapButton = new QPushButton( this->m_ImagesTab );
  copyColormapButton->setText( "Copy Colormap to Other Images" );
  studiesLayout->addWidget( copyColormapButton );

  QObject::connect( copyColormapButton, SIGNAL( clicked() ), this, SLOT( slotCopyColormapToOtherImages() ) );
}

void
QtTriplanarViewer::slotAddStudy( const char* fname )
{
  Study::SmartPtr newStudy( new Study( fname ) );
  this->m_StudiesListBox->insertItem( QString( newStudy->GetFileSystemPath() ) );

  this->m_Studies.push_back( newStudy );
  this->m_ControlsTab->setTabEnabled( this->m_ImagesTab, this->m_Studies.size() > 1 );

  this->slotSwitchToStudy( newStudy );
  this->slotCenter();
}

void
QtTriplanarViewer::slotLoadStudy()
{
  QString path = Q3FileDialog::getExistingDirectory( QString::null, this, "get existing directory", "Load Study" );

  if ( ! (path.isEmpty() || path.isNull() ) ) 
    {
    Study::SmartPtr newStudy( new Study( path.latin1() ) );

    this->m_Studies.push_back( newStudy );
    this->m_ControlsTab->setTabEnabled( this->m_ImagesTab, this->m_Studies.size() > 1 );

    this->m_StudiesListBox->insertItem( QString( newStudy->GetFileSystemPath() ) );
    this->m_StudiesListBox->setCurrentItem( this->m_StudiesListBox->count()-1 );

    this->slotSwitchToStudy( newStudy );
    this->slotCenter();
    }
}

void
QtTriplanarViewer::slotLoadFile()
{
  QString path = Q3FileDialog::getOpenFileName( QString::null, "*", this, "get existing file", "Load File" );
  
  if ( ! (path.isEmpty() || path.isNull() ) ) 
    {
    Study::SmartPtr newStudy( new Study( path.latin1() ) );

    this->m_Studies.push_back( newStudy );
    this->m_ControlsTab->setTabEnabled( this->m_ImagesTab, this->m_Studies.size() > 1 );
    
    this->m_StudiesListBox->insertItem( QString( newStudy->GetFileSystemPath() ) );
    this->m_StudiesListBox->setCurrentItem( this->m_StudiesListBox->count()-1 );

    this->slotSwitchToStudy( newStudy );
    this->slotCenter();
    }
}

void
QtTriplanarViewer::slotReloadData()
{
  if ( this->m_Study ) 
    {
    this->m_Study->ReadVolume( true /* reload */, AnatomicalOrientation::ORIENTATION_STANDARD );
    }
}

void
QtTriplanarViewer::slotSwitchStudy( const QString& study )
{
  for ( size_t idx = 0; idx < this->m_Studies.size(); ++idx )
    {
    if ( study == QString( this->m_Studies[idx]->GetFileSystemPath() ) )
      {
      this->slotSwitchToStudyInternal( this->m_Studies[idx] );
      return;
      }
    }
}

void
QtTriplanarViewer::slotCopyColormapToOtherImages()
{
  if ( this->m_Study ) 
    {
    for ( size_t i = 0; i < this->m_Studies.size(); ++i )
      {
      if ( this->m_Studies[i] != this->m_Study )
	{
	this->m_Studies[i]->CopyColormap( this->m_Study );
	}
      }
    }
}

} // namespace cmtk
