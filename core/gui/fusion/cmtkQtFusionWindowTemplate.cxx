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

#include <cmtkQtFusionWindowTemplate.h>
#include <cmtkQtFusionGlobal.h>

#include <qpushbutton.h>
#include <q3hgroupbox.h>
#include <q3filedialog.h>
#include <Q3VBoxLayout>
#include <Q3PopupMenu>

namespace
cmtk
{

QtFusionWindowTemplate::QtFusionWindowTemplate
( QtSimpleFusionApp *const fusionApp, QWidget *const parent, 
  const char* name, Qt::WFlags flags ) 
  :
  QWidget( parent, name, flags ),
  ZoomFactorPercent( 100 ),
  AxesMode( false ),
  FusionApp( fusionApp )
{
  this->setIcon( QtFusionGlobal::WindowIcon() );

  ViewMenu = new Q3PopupMenu();
  ViewMenu->setCheckable( true );
  ViewMenu->insertItem( "25%", VIEWMENU_25 );
  ViewMenu->insertItem( "33%", VIEWMENU_33 );
  ViewMenu->insertItem( "&50%", VIEWMENU_50 );
  ViewMenu->insertItem( "&66%", VIEWMENU_66 );
  ViewMenu->insertSeparator();
  ViewMenu->insertItem( "&100%", VIEWMENU_100 );
  ViewMenu->insertItem( "&200%", VIEWMENU_200 );
  ViewMenu->insertItem( "&300%", VIEWMENU_300 );
  ViewMenu->insertItem( "&400%", VIEWMENU_400 );
  ViewMenu->setItemChecked( VIEWMENU_100, true );
  ViewMenu->insertSeparator();
  ViewMenu->insertItem( "Annotations...", VIEWMENU_ANNOTATE );
  //  ViewMenu->insertSeparator();
  QObject::connect( ViewMenu, SIGNAL( activated( int ) ), this, SLOT( slotViewMenuCmd( int ) ) );
  
  MenuBar = new QMenuBar( this );
  MenuBar->insertItem( "&View", ViewMenu );
  MenuBar->show();

  MasterLayout = new Q3VBoxLayout( this );
  MasterLayout->setMenuBar( MenuBar );

  WindowLayout = new Q3BoxLayout( MasterLayout, Q3BoxLayout::TopToBottom );

  ViewLayout = new Q3BoxLayout( WindowLayout, Q3BoxLayout::TopToBottom );
  ControlsLayout = new Q3VBoxLayout( WindowLayout );

  ButtonBox = new Q3HGroupBox( this );
  MasterLayout->addWidget( ButtonBox );

  QPushButton* updateButton = new QPushButton( "Update", ButtonBox, "UpdateButton" );
  QObject::connect( updateButton, SIGNAL( clicked() ), SIGNAL( signalUpdate() ) );

  QPushButton* exportButton = new QPushButton( "Export...", ButtonBox, "ExportButton" );

  QPushButton* closeButton = new QPushButton( "Close", ButtonBox, "CloseButton" );

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
QtFusionWindowTemplate::slotViewMenuCmd( int id )
{
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
    qWarning( "Undhandled command in "
	      "QtFusionWindowTemplate::slotViewMenuCmd()." );
    break;
  }

  ViewMenu->setItemChecked( VIEWMENU_25, id == VIEWMENU_25 );
  ViewMenu->setItemChecked( VIEWMENU_33, id == VIEWMENU_33 );
  ViewMenu->setItemChecked( VIEWMENU_50, id == VIEWMENU_50 );
  ViewMenu->setItemChecked( VIEWMENU_66, id == VIEWMENU_66 );
  ViewMenu->setItemChecked( VIEWMENU_100, id == VIEWMENU_100 );
  ViewMenu->setItemChecked( VIEWMENU_200, id == VIEWMENU_200 );
  ViewMenu->setItemChecked( VIEWMENU_300, id == VIEWMENU_300 );
  ViewMenu->setItemChecked( VIEWMENU_400, id == VIEWMENU_400 );

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
  QString path = Q3FileDialog::getSaveFileName( QString::null, "Images (*.pgm *.ppm *.tif)", this, "get destination file", "Save Image" );
  if ( ! (path.isEmpty() || path.isNull() ) )
    this->Export( path );
}

} // namespace cmtk
