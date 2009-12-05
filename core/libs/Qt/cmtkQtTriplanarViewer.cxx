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

#include <qmenu.h>
#include <qmenubar.h>
#include <qapplication.h>
#include <qfiledialog.h>
#include <qlayout.h>
#include <QVBoxLayout>

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
  this->setWindowTitle( "Triplanar Image Viewer" );

  QMenu* StudyMenu = new QMenu();
  StudyMenu->setTitle( "&Study" );
  StudyMenu->addAction( "&Load...", this, SLOT( slotLoadFile() )  );
  StudyMenu->addAction( "&Reload Data...", this, SLOT( slotReloadData() )  );
  StudyMenu->addSeparator();
  StudyMenu->addAction( "&Save" );
  StudyMenu->addAction( "Save &as..." );
  StudyMenu->addAction( "&Export landmarks..." );
  StudyMenu->addSeparator();
  StudyMenu->addAction( "&Quit", qApp, SLOT( quit() ) );
  
  QtImageOperators* imageOperators = new QtImageOperators( &this->m_Study, this, NULL /*progressInstance*/ );
  QObject::connect( imageOperators, SIGNAL( dataChanged( Study::SmartPtr& ) ), this, SLOT( slotDataChanged( Study::SmartPtr& ) ) );
  
  // insert "Study" as first menu.
  MenuBar->insertMenu( this->ViewMenu->menuAction(), StudyMenu );
  // insert "Operators" after "View".
  MenuBar->addMenu( imageOperators->CreatePopupMenu() );
  MenuBar->show();
  
  this->m_ImagesTab = new QWidget( this->m_ControlsTab );
  this->m_ControlsTab->addTab( this->m_ImagesTab, "Images" );
  this->m_ControlsTab->setTabEnabled( this->m_ControlsTab->indexOf( this->m_ImagesTab ), false );

  this->m_StudiesListBox = new QListWidget( this->m_ImagesTab );
  this->m_StudiesListBox->setSelectionMode( QAbstractItemView::SingleSelection );
//  this->m_StudiesListBox->setColumnMode( QListWidget::FixedNumber );
  QObject::connect( this->m_StudiesListBox, SIGNAL( currentTextChanged( const QString& ) ), this, SLOT( slotSwitchStudy( const QString& ) ) );

  QVBoxLayout* studiesLayout = new QVBoxLayout( this->m_ImagesTab );
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
  this->m_StudiesListBox->addItem( QString( newStudy->GetFileSystemPath() ) );

  this->m_Studies.push_back( newStudy );
  this->m_ControlsTab->setTabEnabled( this->m_ControlsTab->indexOf( this->m_ImagesTab ), this->m_Studies.size() > 1 );

  this->slotSwitchToStudy( newStudy );
  this->slotCenter();
}

void
QtTriplanarViewer::slotLoadFile()
{
#ifdef CMTK_BUILD_NRRD
  QString path = QFileDialog::getOpenFileName( this, "Load File", QString(), "All image files (*.hdr *.nii *.nii.gz *.nrrd *.nhdr *.pic);; NIfTI / Analyze (*.hdr *.nii *.nii.gz);; Nrrd (*.nhdr *.nrrd);; BIORAD (*.pic)" );
#else
  QString path = QFileDialog::getOpenFileName( this, "Load File", QString(), "All image files (*.hdr *.nii *.nii.gz *.pic);; NIfTI / Analyze (*.hdr *.nii *.nii.gz);; BIORAD (*.pic)" );
#endif
  
  if ( ! (path.isEmpty() || path.isNull() ) ) 
    {
    Study::SmartPtr newStudy( new Study( path.toLatin1() ) );

    this->m_Studies.push_back( newStudy );
    this->m_ControlsTab->setTabEnabled( this->m_ControlsTab->indexOf( this->m_ImagesTab ), this->m_Studies.size() > 1 );
    
    this->m_StudiesListBox->addItem( QString( newStudy->GetFileSystemPath() ) );
    this->m_StudiesListBox->setCurrentItem( this->m_StudiesListBox->item( this->m_StudiesListBox->count()-1 ) );

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
