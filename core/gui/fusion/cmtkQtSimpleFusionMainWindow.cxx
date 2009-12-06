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
#include <qfiledialog.h>
#include <qinputdialog.h>
#include <qstringlist.h>

#include <iostream>
#include <fstream>

#include <cmtkAffineXform.h>
#include <cmtkClassStreamStudyList.h>
#include <cmtkResourceFile.h>

#include <cmtkQtFusionSlicer.h>
#include <cmtkQtTriplanarViewer.h>
#include <cmtkQtVolumeProperties.h>
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
    RecentListsMenu->addAction( menuItem.sprintf( "&%d. %s", idx, it->c_str() ) );
    ++it;
    ++idx;
    }
}

void
QtSimpleFusionMainWindow::slotRecentListsMenu( QAction* action )
{
  QString path = action->text();
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
  QString path = QFileDialog::getExistingDirectory( this, "Open Studylist" );
  if ( ! (path.isEmpty() || path.isNull() ) )
    this->slotOpenStudyList( path );
}

void 
QtSimpleFusionMainWindow::slotOpenStudyList( const QString& path )
{
  if ( path.isEmpty() || path.isNull() ) 
    {
    QMessageBox::warning( NULL, "Notification", "This is not a valid studylist -\nkeeping previous list instead.", QMessageBox::Ok );
    } 
  else
    {
    StudyList::SmartPtr newStudyList( ClassStreamStudyList::Read( path.toLatin1() ) );
    if ( newStudyList ) 
      {
      FusionApp->slotSetStudyList( newStudyList );
      FusionApp->m_ResourceFile.AddUnique( "RecentLists", path.toLatin1(), 8 );
      this->slotUpdateRecentListsMenu();
      } 
    else
      {
      QMessageBox::critical( NULL, "Error", "Reading of studylist failed.", QMessageBox::Ok );
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
    RecentStudiesMenu->addAction( menuItem.sprintf( "&%d. %s", idx, it->c_str() ) );
    ++it;
    ++idx;
    }
}

void
QtSimpleFusionMainWindow::slotRecentStudiesMenu( QAction* action )
{
  QString path = action->text();
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
  QStringList paths = QFileDialog::getOpenFileNames( this, "Add Study", "*" );
  
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
    QMessageBox::warning( NULL, "Notification", "This is not a valid study.", QMessageBox::Ok );
    } 
  else 
    {
    Study::SmartPtr newStudy( Study::Read( path.toLatin1() ) );
    FusionApp->slotAddStudy( newStudy );
    FusionApp->m_ResourceFile.AddUnique( "RecentStudies", path.toLatin1(), 12 );
    this->slotUpdateRecentStudiesMenu();
    }
}

void
QtSimpleFusionMainWindow::slotSaveStudy()
{
  if ( ! CurrentStudy ) 
    {
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok );
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
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok );
    return;
    }
  
  QString path = QFileDialog::getExistingDirectory( this, "Write Study to" );
  if ( ! (path.isEmpty() || path.isNull() ) ) 
    {
    CurrentStudy->WriteTo( path.toLatin1() );
    }
}

void
QtSimpleFusionMainWindow::slotStudyReadColorMap()
{
  if ( CurrentStudy ) 
    {
    QString path = QFileDialog::getOpenFileName( this, "Read Colormap", QString(), "*.txt" );
    
    if ( ! (path.isEmpty() || path.isNull() ) ) 
      {
      std::ifstream stream( path.toLatin1() );
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
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok );
    return;
    }
  
  QtTriplanarViewer *tpv = new QtTriplanarViewer();
  tpv->slotSwitchToStudy( CurrentStudy );
  QObject::connect( StudyTabs->currentWidget(), SIGNAL( colormap( Study::SmartPtr& ) ), tpv, SLOT( slotColormapChanged( Study::SmartPtr& ) ) );
  
  tpv->show();
}

void 
QtSimpleFusionMainWindow::slotVolumeProperties()
{
  if ( ! CurrentStudy ) 
    {
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok );
    return;
    }
  
  QtVolumeProperties *vp = new QtVolumeProperties( this->CurrentStudy );
  vp->show();
}

void
QtSimpleFusionMainWindow::slotOperatorsMenu( QAction* action )
{
  const int command = action->data().toInt();
  switch ( command ) 
    {
    case OPERATORS_MENU_MEDIAN: 
    {
    if ( CurrentStudy && CurrentStudy->GetVolume() ) 
      {
      bool ok;
      int size = QInputDialog::getInt( this, "Median Filter", "Neighborhood size:",  1, 1, 5, 1, &ok );
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
      int bins = QInputDialog::getInt( this, "Histogram Equalization", "Number of Histogram Bins:",  256, 2, 256, 1, &ok );
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
    if ( CurrentStudy && CurrentStudy->GetVolume() && CurrentStudy->GetVolume()->GetData() ) 
      {
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
QtSimpleFusionMainWindow::slotXformMenuCreate()
{
  if ( ! CurrentStudy ) 
    {
    QMessageBox::warning( NULL, "Notification", "No study currently selected.", QMessageBox::Ok );
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
    QMessageBox::warning( NULL, "Notification", "No target studies are available.", QMessageBox::Ok );
    return;
    }
  
  bool ok = false;
  QString target = QInputDialog::getItem( this, "Select Target Study", "Target", nameList, 0 /*current*/, false /*editable*/, &ok /*ok*/ );
  
  if ( ok ) 
    {
    Study::SmartPtr targetStudy = FusionApp->m_StudyList->FindStudyName( target.toLatin1() );
    
    if ( targetStudy.IsNull() ) 
      {
      QMessageBox::warning( NULL, "Internal Error", "Could not find study with selected name.", QMessageBox::Ok );
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
  
void
QtSimpleFusionMainWindow::slotFusionMenu( QAction* action )
{
  const int command = action->data().toInt();
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
