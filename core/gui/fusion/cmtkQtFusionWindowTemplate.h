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

#ifndef __cmtkQtFusionWindowTemplate_h_included_
#define __cmtkQtFusionWindowTemplate_h_included_

#include <cmtkconfig.h>

#include <qwidget.h>
#include <qlayout.h>
#include <q3groupbox.h>
#include <qmenubar.h>
#include <q3popupmenu.h>

#include <cmtkQtSimpleFusionApp.h>
#include <cmtkStudyList.h>

#include <cmtkQtFusionSlicer.h>
#include <cmtkColormap.h>

namespace
cmtk
{

/** Template class for all fusion windows.
 * This class provides layout boxes for viewers and buttons as well as a basic
 * menu bar.
 */
class QtFusionWindowTemplate :
  /// This is a window widget.
  public QWidget
{
  Q_OBJECT

public:
  /// Constructor.
  QtFusionWindowTemplate( QtSimpleFusionApp *const fusionApp, QWidget *const parent = 0, const char* name = 0, Qt::WFlags flags = 0 );

  /** Destructor.
   * Takes care of destroying the internal visualization pipelines.
   */
  virtual ~QtFusionWindowTemplate();

public slots:
  /// Set fusion slicer object.
  void slotSetFusionSlicer( QtFusionSlicer *const fusionSlicer );

  /// React to changed reference study.
  void slotUpdateReferenceStudy();

  /// Handle "Export" button.
  void slotExport();

  /** Calls the virtual function that writes the current viewer content.
   *\param path File system path for the output file.
   *\param format Image file format.
   @\param comments List of file comments.
  */
  void slotExport( const QString& path, const QString& format, const QStringList* comments )
  {
    this->Export( path, format, comments );
  }

private slots:
  /// React to "View" menu.
  void slotViewMenuCmd( int id );

  /// React to changed slice location etc.
  void slotSliceChanged() { this->UpdateSlice(); }

signals:
  /// Signal emitted by "Update" button.
  void signalUpdate();

  /// Update annotations.
  void updateViewer();

protected:
  /** The current reference study.
   */
  Study::SmartPtr ReferenceStudy;

  /// Zoom factor.
  unsigned int ZoomFactorPercent;

  /// Flag for coordinate axis display.
  bool AxesMode;

  /// The fusion slicer object.
  QtFusionSlicer* FusionSlicer;

  /// Window menu bar.
  QMenuBar* MenuBar;

  /// "View" menu.
  Q3PopupMenu* ViewMenu;

  /// Qt master layout for this window.
  Q3BoxLayout* MasterLayout;

  /// Qt layout for main part of window (excluding buttons).
  Q3BoxLayout* WindowLayout;

  /// Qt layout for top of the window (viewer area).
  Q3BoxLayout* ViewLayout;

  /// Qt layout for center of the window (optional controls area).
  Q3BoxLayout* ControlsLayout;

  /// Qt box for bottom of the window (button area).
  Q3GroupBox* ButtonBox;

  /// This function is called by all slots that change the list of studies.
  virtual void UpdateStudySelection() {};

  /// Export image: This function needs to be overwritten by derived classes.
  virtual void Export( const QString& path, const QString& format = QString::null, const QStringList* comments = NULL ) = 0;

  /// This virtual member is called when the slice changes.
  virtual void UpdateSlice() {};

private:
  /// Command id's for "View" menu.
  typedef enum {
    /// Zoom 25%.
    VIEWMENU_25,
    /// Zoom 33%.
    VIEWMENU_33,
    /// Zoom 50%.
    VIEWMENU_50,
    /// Zoom 66%.
    VIEWMENU_66,
    /// Zoom 100%.
    VIEWMENU_100,
    /// Zoom 200%.
    VIEWMENU_200,
    /// Zoom 300%.
    VIEWMENU_300,
    /// Zoom 400%.
    VIEWMENU_400,
    /// Configure viewer annotations. 
    VIEWMENU_ANNOTATE
  } igsViewMenuCmd;

protected:
  /// Pointer to fusion application.
  QtSimpleFusionApp* FusionApp;
};

} // namespace cmtk

#endif // #ifndef  __cmtkQtFusionWindowTemplate_h_included_

