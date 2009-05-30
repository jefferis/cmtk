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

#include <cmtkQtSimpleFusionMainWindow.h>

#include <cmtkQtSimpleFusionApp.h>

#include <qapplication.h>
#include <qmessagebox.h>
#include <qmenubar.h>
#include <q3textstream.h>
#include <q3filedialog.h>
#include <qinputdialog.h>
#include <qstringlist.h>

#include <iostream>
#include <fstream>

#include <cmtkAffineXform.h>
#include <cmtkClassStreamStudyList.h>
#include <cmtkResourceFile.h>

#include <cmtkQtFusionSlicer.h>
#include <cmtkQtTriplanarViewer.h>
#include <cmtkQtStudyWidget.h>

#include <cmtkQtSeparateView.h>
#include <cmtkQtFusionAlpha.h>
#include <cmtkQtFusionEdge.h>

#include <cmtkSegmentationLabel.h>
#include <cmtkSegmentationLabelIO.h>

namespace
cmtk
{

QtSimpleFusionMainWindow::QtSimpleFusionMainWindow
( QtSimpleFusionApp *const fusionApp )
  : CurrentStudy(),
    QtProgressInstance( NULL )
{
  FusionApp = fusionApp;
  this->InitWidget();
}

QtSimpleFusionMainWindow::~QtSimpleFusionMainWindow()
{
  if ( QtProgressInstance ) delete QtProgressInstance;
}

void
QtSimpleFusionMainWindow::slotUpdateRecentListsMenu()
{
  RecentListsMenu->clear();

  ResourceSection& section = FusionApp->m_ResourceFile["RecentLists"];
  ResourceSection::const_iterator it = section.begin();

  unsigned int idx = 1;
  while ( it != section.end() ) 
    {
    QString menuItem;
    RecentListsMenu->insertItem( menuItem.sprintf( "&%d. %s", idx, it->c_str() ) );
    ++it;
    ++idx;
    }
}

void
QtSimpleFusionMainWindow::slotRecentListsMenu( const int id )
{
  QString path = RecentListsMenu->text( id );
  if ( path.length() > 4 ) 
    {
    // open studylist referenced by Recent menu entry; chop of first
    // four characters ("&n. " prefix)
    this->slotOpenStudyList( path.right( path.length()-4 ) );
    }
}

void 
QtSimpleFusionMainWindow::slotOpenStudyList()
{
  QString path = Q3FileDialog::getExistingDirectory( QString::null, this, "get existing directory", "Open Studylist" );
  if ( ! (path.isEmpty() || path.isNull() ) )
    this->slotOpenStudyList( path );
}

void 
QtSimpleFusionMainWindow::slotOpenStudyList( const QString& path )
{
  if ( path.isEmpty() || path.isNull() ) 
    {
    QMessageBox::warning( NULL, "Notification", "This is not a valid studylist -\nkeeping previous list instead.", 
			  QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
    } 
  else
    {
    StudyList::SmartPtr newStudyList( ClassStreamStudyList::Read( path.latin1() ) );
    if ( newStudyList ) 
      {
      FusionApp->slotSetStudyList( newStudyList );
      FusionApp->m_ResourceFile.AddUnique( "RecentLists", path.latin1(), 8 );
      this->slotUpdateRecentListsMenu();
      } 
    else
      {
      QMessageBox::critical( NULL, "Error", "Reading of studylist failed.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
      }
    }
}

void
QtSimpleFusionMainWindow::slotUpdateRecentStudiesMenu()
{
  RecentStudiesMenu->clear();

  ResourceSection& section = FusionApp->m_ResourceFile["RecentStudies"];
  ResourceSection::const_iterator it = section.begin();

  unsigned int idx = 1;
  while ( it != section.end() ) 
    {
    QString menuItem;
    RecentStudiesMenu->insertItem( menuItem.sprintf( "&%d. %s", idx, it->c_str() ) );
    ++it;
    ++idx;
    }
}

void
QtSimpleFusionMainWindow::slotRecentStudiesMenu( const int id )
{
  QString path = RecentStudiesMenu->text( id );
  if ( path.length() > 4 ) 
    {
    // open studylist referenced by Recent menu entry; chop of first
    // four characters ("&n. " prefix)
    this->slotAddStudy( path.right( path.length()-4 ) );
    }
}

void 
QtSimpleFusionMainWindow::slotAddStudy()
{
  QString path = Q3FileDialog::getExistingDirectory( QString::null, this, "get existing directory", "Add Study" );
  if ( ! (path.isEmpty() || path.isNull() ) )
    this->slotAddStudy( path );
}

void 
QtSimpleFusionMainWindow::slotAddStudyFiles()
{
  QStringList paths = Q3FileDialog::getOpenFileNames( QString::null, "*", this, "get existing file", "Add Study" );
  
  QStringList::const_iterator it = paths.begin();
  while ( it != paths.end() ) 
    {
    if ( ! ((*it).isEmpty() || (*it).isNull() ) )
      this->slotAddStudy( *it );
    ++it;
    }
}

void 
QtSimpleFusionMainWindow::slotAddStudy( const QString& path )
{
  if ( path.isEmpty() || path.isNull() ) 
    {
    QMessageBox::warning( NULL, "Notification", "This is not a valid study.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
    } 
  else 
    {
    Study::SmartPtr newStudy( Study::Read( path.latin1() ) );
    FusionApp->slotAddStudy( newStudy );
    FusionApp->m_ResourceFile.AddUnique( "RecentStudies", path.latin1(), 12 );
    this->slotUpdateRecentStudiesMenu();
    }
}

void
QtSimpleFusionMainWindow::slotSaveStudy()
{
  if ( ! CurrentStudy ) 
    {
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
    return;
    }
  
  if ( ! CurrentStudy->GetFileSystemPath() )
    this->slotSaveStudyAs();
  else
    CurrentStudy->Write();
}

void
QtSimpleFusionMainWindow::slotSaveStudyAs()
{
  if ( ! CurrentStudy ) 
    {
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
    return;
    }
  
  QString path = Q3FileDialog::getExistingDirectory( QString::null, this, "get existing directory", "Write Study to" );
  if ( ! (path.isEmpty() || path.isNull() ) ) 
    {
    CurrentStudy->WriteTo( path.latin1() );
    }
}

void
QtSimpleFusionMainWindow::slotStudyReadColorMap()
{
  if ( CurrentStudy ) 
    {
    QString path = Q3FileDialog::getOpenFileName( QString::null, "*.txt", this, "get existing file", "Read Colormap" );
    
    if ( ! (path.isEmpty() || path.isNull() ) ) 
      {
      std::ifstream stream( path.latin1() );
      SegmentationLabelMap userColorMap;
      stream >> userColorMap;
      CurrentStudy->SetFromLabelMap( userColorMap );
      }
    }
}

void
QtSimpleFusionMainWindow::slotStudyReload()
{
  if ( CurrentStudy ) 
    {
    CurrentStudy->ReadVolume( true /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
    }
}

void 
QtSimpleFusionMainWindow::slotTriplanarViewer()
{
  if ( ! CurrentStudy ) 
    {
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
    return;
    }
  
  QtTriplanarViewer *tpv = new QtTriplanarViewer();
  tpv->slotSwitchToStudy( CurrentStudy );
  QObject::connect( StudyTabs->currentPage(), SIGNAL( colormap( Study::SmartPtr& ) ), tpv, SLOT( slotColormapChanged( Study::SmartPtr& ) ) );
  
  tpv->show();
}

void
QtSimpleFusionMainWindow::slotOperatorsMenu( int command )
{
  switch ( command ) 
    {
    case OPERATORS_MENU_MEDIAN: 
    {
    if ( CurrentStudy && CurrentStudy->GetVolume() ) 
      {
      bool ok;
      int size = QInputDialog::getInteger( "Median Filter", "Neighborhood size:",  1, 1, 5, 1, &ok, this );
      if ( ok ) 
	{
        // user entered something and pressed OK
	QtProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_DIALOG );
	CurrentStudy->GetVolume()->ApplyMedianFilter( 1 + 2 * size );
	CurrentStudy->UpdateFromVolume();
	FusionApp->slotDataChanged( CurrentStudy );
	} 
      else 
	{
        // user pressed Cancel
	}
      }
    break;
    }
    case OPERATORS_MENU_SOBEL: 
    {
    if ( CurrentStudy && CurrentStudy->GetVolume() ) 
      {
      QtProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_BAR );
      CurrentStudy->GetVolume()->ApplySobelFilter();
      CurrentStudy->UpdateFromVolume();
      }
    break;
    }
    case OPERATORS_MENU_HISTOGRAM: 
    {
    if ( CurrentStudy && CurrentStudy->GetVolume() && CurrentStudy->GetVolume()->GetData() ) 
      {
      QtProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_BAR );
      bool ok;
      int bins = QInputDialog::getInteger( "Histogram Equalization", "Number of Histogram Bins:",  256, 2, 256, 1, &ok, this );
      if ( ok )
	{
	// user entered something and pressed OK
	QtProgressInstance->SetProgressWidgetMode( QtProgress::PROGRESS_DIALOG );
	CurrentStudy->GetVolume()->GetData()->HistogramEqualization( bins );
	CurrentStudy->UpdateFromVolume();
	FusionApp->slotDataChanged( CurrentStudy );
	} 
      else 
	{
        // user pressed Cancel
	}
      }
    break;
    case OPERATORS_MENU_ABS: 
    case OPERATORS_MENU_LOG: 
    case OPERATORS_MENU_EXP: 
    {
    if ( CurrentStudy && CurrentStudy->GetVolume() && CurrentStudy->GetVolume()->GetData() ) {
    switch ( command )
      {
      case OPERATORS_MENU_ABS: 
	CurrentStudy->GetVolume()->GetData()->ApplyFunction( fabs );
	CurrentStudy->UpdateFromVolume();
	FusionApp->slotDataChanged( CurrentStudy );
	break;
      case OPERATORS_MENU_LOG: 
	CurrentStudy->GetVolume()->GetData()->ApplyFunction( log );
	CurrentStudy->UpdateFromVolume();
	FusionApp->slotDataChanged( CurrentStudy );
	break;
      case OPERATORS_MENU_EXP: 
	CurrentStudy->GetVolume()->GetData()->ApplyFunction( exp );
	CurrentStudy->UpdateFromVolume();
	FusionApp->slotDataChanged( CurrentStudy );
	break;
      }
    }
    break;
    }
    }
    }
}

void
QtSimpleFusionMainWindow::slotXformMenu( int command )
{
  switch ( command ) 
    {
    case MENU_XFORM_CREATE: 
    {
    if ( ! CurrentStudy ) 
      {
      QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
      return;
      }
    
    std::list<Study::SmartPtr> targets = FusionApp->m_StudyList->GetIndependentStudyList( CurrentStudy );
    
    QStringList nameList;
    std::list<Study::SmartPtr>::const_iterator it = targets.begin();
    while ( it != targets.end() ) 
      {
      nameList.push_back( (*it)->GetName() );
      ++it;
      }
    
    if ( ! nameList.size() ) 
      {
      QMessageBox::warning( NULL, "Notification", "No target studies are available.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
      return;
      }
    
    bool ok = false;
    QString target = QInputDialog::getItem( "Select Target Study", "Target", nameList, 0 /*current*/, false /*editable*/, &ok /*ok*/, this );
    
    if ( ok ) 
      {
      Study::SmartPtr targetStudy = FusionApp->m_StudyList->FindStudyName( target.latin1() );
      
      if ( targetStudy.IsNull() ) 
	{
	QMessageBox::warning( NULL, "Internal Error", "Could not find study with selected name.", 
			      QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
	return;
	}
      
      AffineXform::SmartPtr newXform( new AffineXform );
      AffineXform::SmartPtr inverse = newXform->GetInverse();
      WarpXform::SmartPtr nullWarp( NULL );
      
      // add identity transformation
      FusionApp->m_StudyList->AddXform( CurrentStudy, targetStudy, newXform, nullWarp );
      // add inverse in opposite direction
      FusionApp->m_StudyList->AddXform( targetStudy, CurrentStudy, inverse, nullWarp );
      
      // notify application of change to studylist
      FusionApp->slotSetStudyList( FusionApp->m_StudyList );
      }
    }
    }
}

void
QtSimpleFusionMainWindow::slotFusionMenu( int command )
{
  switch ( command ) 
    {
    case FUSION_MENU_SEPARATE: 
    {
    QtSeparateView* separateView = new QtSeparateView( FusionApp );
    QObject::connect( this->FusionApp->GetFusionSlicer(), SIGNAL( signalExport( const QString&, const QString&, const QStringList* ) ),
		      separateView, SLOT( slotExport( const QString&, const QString&, const QStringList* ) ) );
    break;
    }
    case FUSION_MENU_OVERLAY:
      break;
    case FUSION_MENU_ALPHA: 
    {
    QtFusionAlpha* fusionAlpha = new QtFusionAlpha( FusionApp );
    QObject::connect( this->FusionApp->GetFusionSlicer(), SIGNAL( signalExport( const QString&, const QString&, const QStringList* ) ),
		      fusionAlpha, SLOT( slotExport( const QString&, const QString&, const QStringList* ) ) );
    break;
    }
    break;
    case FUSION_MENU_EDGE: 
    {
    QtFusionEdge* fusionEdge = new QtFusionEdge( FusionApp );
    QObject::connect( this->FusionApp->GetFusionSlicer(), SIGNAL( signalExport( const QString&, const QString&, const QStringList* ) ),
		      fusionEdge, SLOT( slotExport( const QString&, const QString&, const QStringList* ) ) );
    break;
    }
    case FUSION_MENU_ISOLINES:
      break;
    case FUSION_MENU_DIFFERENCE:
      break;
    case FUSION_MENU_ROI:
      break;
    case FUSION_MENU_MIX:
      break;
    case FUSION_MENU_COLOR:
      break;
    case FUSION_MENU_SLICER: 
    {
    QtFusionSlicer* fusionSlicer = FusionApp->GetFusionSlicer();
    fusionSlicer->show();
    break;
    }
    default:
      qWarning( "In QtSimpleFusionMainWindow::slotOpenFusionWindow: command %d not handled.", command );
      break;
    }
}

void
QtSimpleFusionMainWindow::slotSwitchStudy( QWidget* currentTab )
{
  QtStudyWidget* studyWidget =dynamic_cast<QtStudyWidget*>( currentTab );
  
  if ( ! studyWidget ) 
    {
    qWarning( "In QtSimpleFusionMainWindow::slotSwitchStudyTab: RTTI says current tab is not QtStudyWidget" );
    CurrentStudy = Study::SmartPtr::Null;
    return;
    }
  
  CurrentStudy = studyWidget->GetStudy();
}

} // namespace cmtk
