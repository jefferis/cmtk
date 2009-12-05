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

#include <cmtkQtFusionWindowTemplate.h>
#include <cmtkQtFusionGlobal.h>

#include <qpushbutton.h>
#include <qgroupbox.h>
#include <qfiledialog.h>
#include <QVBoxLayout>
#include <QMenu>

namespace
cmtk
{

QtFusionWindowTemplate::QtFusionWindowTemplate
( QtSimpleFusionApp *const fusionApp, QWidget *const parent, Qt::WFlags flags ) 
  :
  QWidget( parent, flags ),
  ZoomFactorPercent( 100 ),
  AxesMode( false ),
  FusionApp( fusionApp )
{
  this->setWindowIcon( QtFusionGlobal::WindowIcon() );

  MenuBar = new QMenuBar( this );
  MenuBar->show();

  ViewMenu = MenuBar->addMenu( "&View" );
  ViewMenu->addAction( "25%" )->setData( QVariant( VIEWMENU_25 ) );
  ViewMenu->addAction( "33%" )->setData( QVariant( VIEWMENU_33 ) );
  ViewMenu->addAction( "&50%" )->setData( QVariant( VIEWMENU_50 ) );
  ViewMenu->addAction( "&66%" )->setData( QVariant( VIEWMENU_66 ) );
  ViewMenu->addSeparator();
  ViewMenu->addAction( "&100%" )->setData( QVariant( VIEWMENU_100 ) );
  ViewMenu->addAction( "&200%" )->setData( QVariant( VIEWMENU_200 ) );
  ViewMenu->addAction( "&300%" )->setData( QVariant( VIEWMENU_300 ) );
  ViewMenu->addAction( "&400%" )->setData( QVariant( VIEWMENU_400 ) );
  ViewMenu->addSeparator();
  ViewMenu->addAction( "Annotations..." )->setData( QVariant( VIEWMENU_ANNOTATE ) );
  QObject::connect( ViewMenu, SIGNAL( activated( QAction* ) ), this, SLOT( slotViewMenuCmd( QAction* ) ) );
  
  MasterLayout = new QVBoxLayout( this );
  MasterLayout->setMenuBar( MenuBar );

  WindowLayout = new QBoxLayout( MasterLayout, QBoxLayout::TopToBottom );

  ViewLayout = new QBoxLayout( WindowLayout, QBoxLayout::TopToBottom );
  ControlsLayout = new QVBoxLayout( WindowLayout );

  ButtonBox = new QGroupBox( this );
  MasterLayout->addWidget( ButtonBox );

  QPushButton* updateButton = new QPushButton( "Update", ButtonBox );
  QObject::connect( updateButton, SIGNAL( clicked() ), SIGNAL( signalUpdate() ) );

  QPushButton* exportButton = new QPushButton( "Export...", ButtonBox );

  QPushButton* closeButton = new QPushButton( "Close", ButtonBox );

  QObject::connect( exportButton, SIGNAL( clicked() ), this, SLOT( slotExport() ) );
  QObject::connect( closeButton, SIGNAL( clicked() ), this, SLOT( close() ) );

  // link to other parts of the application
  FusionSlicer = FusionApp->GetFusionSlicer();
  FusionSlicer->show();

  QObject::connect( FusionSlicer, SIGNAL( sliceChanged() ), this, SLOT( slotSliceChanged() ) );
  QObject::connect( FusionApp, SIGNAL( signalReferenceStudyChanged() ), this, SLOT( slotUpdateReferenceStudy() ) );
}

QtFusionWindowTemplate::~QtFusionWindowTemplate()
{
}

void 
QtFusionWindowTemplate::slotViewMenuCmd( QAction* action )
{
  const int id = action->data().toInt();
  switch ( id ) {
  case VIEWMENU_25:
    ZoomFactorPercent = 25;
    break;
  case VIEWMENU_33:
    ZoomFactorPercent = 33;
    break;
  case VIEWMENU_50:
    ZoomFactorPercent = 50;
    break;
  case VIEWMENU_66:
    ZoomFactorPercent = 66;
    break;
  case VIEWMENU_100:
    ZoomFactorPercent = 100;
    break;
  case VIEWMENU_200:
    ZoomFactorPercent = 200;
    break;
  case VIEWMENU_300:
    ZoomFactorPercent = 300;
    break;
  case VIEWMENU_400:
    ZoomFactorPercent = 400;
    break;
  case VIEWMENU_ANNOTATE:
    break;
  default:
    qWarning( "Undhandled command in QtFusionWindowTemplate::slotViewMenuCmd()." );
    break;
  }

  emit updateViewer();
}

void 
QtFusionWindowTemplate::slotSetFusionSlicer
( QtFusionSlicer *const fusionSlicer )
{
  FusionSlicer = fusionSlicer;
  this->UpdateStudySelection();
}

void
QtFusionWindowTemplate::slotUpdateReferenceStudy()
{
  this->UpdateStudySelection();
}

void
QtFusionWindowTemplate::slotExport()
{
  QString path = QFileDialog::getSaveFileName( QString::null, "Images (*.pgm *.ppm *.tif)", this, "get destination file", "Save Image" );
  if ( ! (path.isEmpty() || path.isNull() ) )
    this->Export( path );
}

} // namespace cmtk
