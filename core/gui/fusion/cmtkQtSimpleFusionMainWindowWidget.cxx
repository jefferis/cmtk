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

#include <cmtkconfig.h>

#include "cmtkQtSimpleFusionMainWindow.h"
#include "cmtkQtFusionGlobal.h"

#include <qapplication.h>
#include <qmessagebox.h>
#include <qmenubar.h>
#include <qstatusbar.h>
#include <qprogressbar.h>
#include <QLabel>
#include <QMenu>

namespace
cmtk
{

void
QtSimpleFusionMainWindow::InitWidget()
{
  this->setWindowTitle( "CMTK Image Fusion" );
  this->setWindowIcon( QtFusionGlobal::WindowIcon() );

  QStatusBar* statBar = this->statusBar();

  StatusLabel = new QLabel( statBar );
  StatusLabel->setText( "Image Fusion" );
  statBar->addWidget( StatusLabel, 4 );

  QProgressBar* progressBar = new QProgressBar( statBar );
  statBar->addWidget( progressBar, 1 );
  progressBar->show();

  QtProgressInstance = new QtProgress( this );
  QtProgressInstance->SetProgressBar( progressBar );
  Progress::SetProgressInstance( QtProgressInstance );

  MenuBar = this->menuBar();

  ListMenu = MenuBar->addMenu( "&List" );
  ListMenu->addAction( "&Open...", this, SLOT( slotOpenStudyList() ) );
  ListMenu->addSeparator();
  ListMenu->addAction( "&Save", this, SLOT( slotOpenStudyList() ) );
  ListMenu->addAction( "Save &As...", this, SLOT( slotOpenStudyList() ) );
  ListMenu->addSeparator();

  RecentListsMenu = ListMenu->addMenu( "&Recent" );
  ListMenu->addSeparator();
  this->slotUpdateRecentListsMenu();
  QObject::connect( RecentListsMenu, SIGNAL( triggered( QAction* ) ), this, SLOT( slotRecentListsMenu( QAction* ) ) );

  ListMenu->addAction( "&Quit", qApp, SLOT( quit() ) );

  StudyMenu = MenuBar->addMenu( "&Study" );
  StudyMenu->addAction( "&Add...", this, SLOT( slotAddStudy() ) );

  RecentStudiesMenu = StudyMenu->addMenu( "&Recent" );
  this->slotUpdateRecentStudiesMenu();
  QObject::connect( RecentStudiesMenu, SIGNAL( triggered( QAction* ) ), this, SLOT( slotRecentStudiesMenu( QAction* ) ) );

  StudyMenu->addSeparator();
  StudyMenu->addAction( "Load &Colormap...", this, SLOT( slotStudyReadColorMap() ) );
  StudyMenu->addAction( "Re&load", this, SLOT( slotStudyReload() ) );
  StudyMenu->addAction( "&Properties...", this, SLOT( slotVolumeProperties() ) );
  StudyMenu->addAction( "&Triplanar Viewer...", this, SLOT( slotTriplanarViewer() ) );

  OperatorsMenu = MenuBar->addMenu( "&Operators" );
  OperatorsMenu->addAction( "&Median Filter..." )->setData( OPERATORS_MENU_MEDIAN );
  OperatorsMenu->addAction( "&Histogram Equalization..." )->setData( OPERATORS_MENU_HISTOGRAM );
  OperatorsMenu->addAction( "&Sobel Edge Filter" )->setData( OPERATORS_MENU_SOBEL );
  OperatorsMenu->addSeparator();
  QObject::connect( OperatorsMenu, SIGNAL( triggered( QAction* ) ), this, SLOT( slotOperatorsMenu( QAction* ) ) );

  AlgOperatorsMenu = OperatorsMenu->addMenu( "&Algebraic" );
  AlgOperatorsMenu->addAction( "&abs()" )->setData( OPERATORS_MENU_ABS );
  AlgOperatorsMenu->addAction( "&log()" )->setData( OPERATORS_MENU_LOG );
  AlgOperatorsMenu->addAction( "&exp()" )->setData( OPERATORS_MENU_EXP );
  QObject::connect( AlgOperatorsMenu, SIGNAL( triggered( QAction* ) ), this, SLOT( slotOperatorsMenu( QAction* ) ) );
  
  XformMenu = MenuBar->addMenu( "&Transform" );;
  XformMenu->addAction( "&Create...", this, SLOT( slotXformMenuCreate() ) );

  FusionMenu = MenuBar->addMenu( "&Fusion" );
  FusionMenu->addAction( "&Separate View" )->setData(  FUSION_MENU_SEPARATE );
  FusionMenu->addAction( "&Overlay" )->setData( FUSION_MENU_OVERLAY );
  FusionMenu->addAction( "&Alpha Blending" )->setData( FUSION_MENU_ALPHA );
  FusionMenu->addAction( "&Edge Blending" )->setData( FUSION_MENU_EDGE );
  FusionMenu->addAction( "&Isolines" )->setData( FUSION_MENU_ISOLINES );
  FusionMenu->addAction( "&Difference" )->setData( FUSION_MENU_DIFFERENCE );
  FusionMenu->addAction( "&Region-of-Interest" )->setData( FUSION_MENU_ROI );
  FusionMenu->addAction( "&Mix" )->setData( FUSION_MENU_MIX );
  FusionMenu->addAction( "&Color/Brightness" )->setData( FUSION_MENU_COLOR );
  FusionMenu->addSeparator();
  FusionMenu->addAction( "&Planar Slicer..." )->setData( FUSION_MENU_SLICER );
  QObject::connect( FusionMenu, SIGNAL( triggered( QAction* ) ), this, SLOT( slotFusionMenu( QAction* ) ) );

  StudyTabs = new QTabWidget( this );
  QObject::connect( StudyTabs, SIGNAL( currentChanged ( QWidget* ) ), this, SLOT( slotSwitchStudy( QWidget* ) ) );
  StudyTabs->show();

  this->setCentralWidget( StudyTabs );
}

} // namespace cmtk
