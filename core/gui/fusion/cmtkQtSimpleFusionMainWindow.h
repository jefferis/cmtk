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

#ifndef __cmtkQtSimpleFusionMainWindow_h_included_
#define __cmtkQtSimpleFusionMainWindow_h_included_

#include <cmtkconfig.h>

#include <qwidget.h>
#include <qlabel.h>
#include <qobject.h>
#include <qmainwindow.h>
#include <qpoint.h>
#include <qtabwidget.h>

#include <cmtkStudy.h>
#include <cmtkStudyList.h>

#include <cmtkQtProgress.h>
#include <cmtkQtFusionSlicer.h>
#include <QMenu>

namespace
cmtk
{

/// Forward declaration.
class QtSimpleFusionApp;

/** Main window class for the Qt image fusion application.
 */
class QtSimpleFusionMainWindow :
  /// Inherit from Qt's MainWindow class.
  public QMainWindow 
{
  Q_OBJECT // we're using slots

signals:
    /// This signal is emitted whenever the colormap changes.
  void colormapChanged( Study::SmartPtr& study );

private slots:
  /// Open studylist with dialog.
  void slotOpenStudyList();

  /// Open given studylist.
  void slotOpenStudyList( const QString& path );

  /// Slot to perform operations from "Transform" menu.
  void slotXformMenuCreate();

  /// Slot to open selected type of fusion window.
  void slotFusionMenu( QAction* );

  /// Slot to open selected type of fusion window.
  void slotOperatorsMenu( QAction* );

  /// Update menu with recent studylists.
  void slotUpdateRecentListsMenu();

  /// This lot handles clicks from the "Recent" studylists menu.
  void slotRecentListsMenu( QAction *action );

  /// Add study with dialog.
  void slotAddStudy();

  /// Add study with dialog.
  void slotAddStudyFiles();

  /// Add given study.
  void slotAddStudy( const QString& path );

  /// Save study to its current location.
  void slotSaveStudy();

  /// Save study to new location.
  void slotSaveStudyAs();

  /// Reload study data.
  void slotStudyReadColorMap();

  /// Reload study data.
  void slotStudyReload();

  /// Add study with dialog.
  void slotTriplanarViewer();

  /// Add study with dialog.
  void slotVolumeProperties();

  /// Update menu with recent studies.
  void slotUpdateRecentStudiesMenu();

  /// This lot handles clicks from the "Recent" studies menu.
  void slotRecentStudiesMenu( QAction* action );

public slots:
  /// Slot to set new study.
  void slotAddStudy( Study::SmartPtr& study );

  /// Slot to delete one study.
  void slotDeleteStudy( Study::SmartPtr& study );

  /// Slot to delete all studies.
  void slotDeleteAllStudies();

  /// Slot called when switching to different study in main widget.
  void slotSwitchStudy( QWidget* currentTab );

public:
  /// Constructor.
  QtSimpleFusionMainWindow( QtSimpleFusionApp *const fusionApp );

  /// Destructor.
  virtual ~QtSimpleFusionMainWindow();

private:
  /// The study currently displayed.
  Study::SmartPtr CurrentStudy;

  /// Initialize the widget components of this object.
  void InitWidget();

  /// The application object.
  QtSimpleFusionApp* FusionApp;

  /// The Qt progress instance.
  QtProgress* QtProgressInstance;

  /// The status bar text label.
  QLabel* StatusLabel;

  /// Tab widget for the studies.
  QTabWidget* StudyTabs;

  /// The menu bar.
  QMenuBar* MenuBar;

  /// The "List" menu.
  QMenu* ListMenu;

  /// Sub-menu with recently opened lists.
  QMenu* RecentListsMenu;

  /// The "Study" menu.
  QMenu* StudyMenu;

  /// Sub-menu with recently opened studies.
  QMenu* RecentStudiesMenu;

  /// The "Operators" menu.
  QMenu* OperatorsMenu;

  /// The "Operators->Algebraic" submenu.
  QMenu* AlgOperatorsMenu;

  /// The "Transformation" menu.
  QMenu* XformMenu;

  /// Enum for commands in the "Study" menu.
  enum {
    REGISTRATION_MENU_MANUAL,
    REGISTRATION_MENU_RIGID
  };

  /// Enum for commands in the "Transform" menu.
  enum {
    MENU_XFORM_CREATE
  };

  /// Enum for commands in the "Operators" menu.
  enum {
    OPERATORS_MENU_MEDIAN,
    OPERATORS_MENU_SOBEL,
    OPERATORS_MENU_HISTOGRAM,
    OPERATORS_MENU_ABS,
    OPERATORS_MENU_LOG,
    OPERATORS_MENU_EXP
  };

  /// Enum for commands in the "Fusion" menu.
  enum {
    FUSION_MENU_SEPARATE,
    FUSION_MENU_OVERLAY,
    FUSION_MENU_ALPHA,
    FUSION_MENU_EDGE,
    FUSION_MENU_ISOLINES,
    FUSION_MENU_DIFFERENCE,
    FUSION_MENU_ROI,
    FUSION_MENU_MIX,
    FUSION_MENU_COLOR,
    FUSION_MENU_3D,
    FUSION_MENU_SLICER
  };

  /// The "Fusion" menu.
  QMenu* FusionMenu;
};

} // namespace cmtk

#endif // #ifndef __cmtkQtSimpleFusionMainWindow_h_included_
