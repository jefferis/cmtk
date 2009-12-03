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

#include <cmtkconfig.h>

#include <cmtkQtSimpleFusionMainWindow.h>
#include <cmtkQtFusionGlobal.h>

#include <qapplication.h>
#include <qmessagebox.h>
#include <qmenubar.h>
#include <q3textstream.h>
#include <q3filedialog.h>
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
  this->setCaption( "IGL Image Fusion" );
  this->setIcon( QtFusionGlobal::WindowIcon() );

  QStatusBar* statBar = this->statusBar();

  StatusLabel = new QLabel( statBar, "StatusLabel" );
  StatusLabel->setText( "Image Fusion" );
  statBar->addWidget( StatusLabel, 4 );

  QProgressBar* progressBar = new QProgressBar( statBar );
  statBar->addWidget( progressBar, 1 );
  progressBar->show();

  QtProgressInstance = new QtProgress( this );
  QtProgressInstance->SetProgressBar( progressBar );
  Progress::SetProgressInstance( QtProgressInstance );

  MenuBar = this->menuBar();

  ListMenu = new QMenu;
  ListMenu->insertItem( "&Open...", this, SLOT( slotOpenStudyList() ) );
  ListMenu->insertSeparator();
  ListMenu->insertItem( "&Save", this, SLOT( slotOpenStudyList() ) );
  ListMenu->insertItem( "Save &As...", this, SLOT( slotOpenStudyList() ) );
  ListMenu->insertSeparator();

  RecentListsMenu = new QMenu;
  ListMenu->insertItem( "&Recent", RecentListsMenu );
  ListMenu->insertSeparator();
  this->slotUpdateRecentListsMenu();
  QObject::connect( RecentListsMenu, SIGNAL( activated( int ) ), this, SLOT( slotRecentListsMenu( int ) ) );

  ListMenu->insertItem( "&Quit", qApp, SLOT( quit() ) );

  RecentStudiesMenu = new QMenu;
  this->slotUpdateRecentStudiesMenu();
  QObject::connect( RecentStudiesMenu, SIGNAL( activated( int ) ), this, SLOT( slotRecentStudiesMenu( int ) ) );

  StudyMenu = new QMenu;
  StudyMenu->insertItem( "&Add...", this, SLOT( slotAddStudy() ) );
  StudyMenu->insertItem( "Add &File...", this, SLOT( slotAddStudyFiles() ) );
  StudyMenu->insertItem( "&Recent", RecentStudiesMenu );
  StudyMenu->insertSeparator();
  StudyMenu->insertItem( "&Save", this, SLOT( slotSaveStudy() ) );
  StudyMenu->insertItem( "Save as...", this, SLOT( slotSaveStudyAs() ) );
  StudyMenu->insertSeparator();
  StudyMenu->insertItem( "Load &Colormap...", this, SLOT( slotStudyReadColorMap() ) );
  StudyMenu->insertItem( "&Edit...", 0 );
  StudyMenu->insertItem( "Re&load", this, SLOT( slotStudyReload() ) );
  StudyMenu->insertItem( "&Delete", 0 );
  StudyMenu->insertItem( "&Properties...", 0 );
  StudyMenu->insertItem( "&Triplanar Viewer...", this, SLOT( slotTriplanarViewer() ) );

  AlgOperatorsMenu = new QMenu;
  AlgOperatorsMenu->insertItem( "&abs()", OPERATORS_MENU_ABS );
  AlgOperatorsMenu->insertItem( "&log()", OPERATORS_MENU_LOG );
  AlgOperatorsMenu->insertItem( "&exp()", OPERATORS_MENU_EXP );
  QObject::connect( AlgOperatorsMenu, SIGNAL( activated( int ) ), this, SLOT( slotOperatorsMenu( int ) ) );

  OperatorsMenu = new QMenu;
  OperatorsMenu->insertItem( "&Median Filter...", OPERATORS_MENU_MEDIAN );
  OperatorsMenu->insertItem( "&Histogram Equalization...", OPERATORS_MENU_HISTOGRAM );
  OperatorsMenu->insertItem( "&Sobel Edge Filter", OPERATORS_MENU_SOBEL );
  OperatorsMenu->insertSeparator();
  OperatorsMenu->insertItem( "&Algebraic", AlgOperatorsMenu );
  QObject::connect( OperatorsMenu, SIGNAL( activated( int ) ), this, SLOT( slotOperatorsMenu( int ) ) );

  
  XformMenu = new QMenu;
  XformMenu->insertItem( "&Create...", MENU_XFORM_CREATE );
  XformMenu->insertSeparator();
  XformMenu->insertItem( "&Affine Transformation Editor...", 0 );
  XformMenu->insertItem( "&Rigid/Affine Registration...", 0 );
  XformMenu->insertItem( "&Visualize Deformation...", 0 );
  QObject::connect( XformMenu, SIGNAL( activated( int ) ), this, SLOT( slotXformMenu( int ) ) );

  FusionMenu = new QMenu;
  FusionMenu->insertItem( "&Separate View", FUSION_MENU_SEPARATE );
  FusionMenu->insertItem( "&Overlay", FUSION_MENU_OVERLAY );
  FusionMenu->insertItem( "&Alpha Blending", FUSION_MENU_ALPHA );
  FusionMenu->insertItem( "&Edge Blending", FUSION_MENU_EDGE );
  FusionMenu->insertItem( "&Isolines", FUSION_MENU_ISOLINES );
  FusionMenu->insertItem( "&Difference", FUSION_MENU_DIFFERENCE );
  FusionMenu->insertItem( "&Region-of-Interest", FUSION_MENU_ROI );
  FusionMenu->insertItem( "&Mix", FUSION_MENU_MIX );
  FusionMenu->insertItem( "&Color/Brightness", FUSION_MENU_COLOR );
  FusionMenu->insertSeparator();
#ifdef IGS_HAVE_VTK
  FusionMenu->insertItem( "&3D View", FUSION_MENU_3D );
  FusionMenu->insertSeparator();
#endif // #ifdef IGS_HAVE_VTK
  FusionMenu->insertItem( "&Planar Slicer...", FUSION_MENU_SLICER );
  QObject::connect( FusionMenu, SIGNAL( activated( int ) ), this, SLOT( slotFusionMenu( int ) ) );

  MenuBar->insertItem( "&List", ListMenu );
  MenuBar->insertItem( "&Study", StudyMenu );
  MenuBar->insertItem( "&Operators", OperatorsMenu );
  MenuBar->insertItem( "&Transform", XformMenu );
  MenuBar->insertItem( "&Fusion", FusionMenu );

  StudyTabs = new QTabWidget( this );
  QObject::connect( StudyTabs, SIGNAL( currentChanged ( QWidget* ) ), this, SLOT( slotSwitchStudy( QWidget* ) ) );
  StudyTabs->show();

  this->setCentralWidget( StudyTabs );
}

} // namespace cmtk
