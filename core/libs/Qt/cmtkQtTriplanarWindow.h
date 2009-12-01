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

#ifndef __cmtkQtTriplanarWindow_h_included_
#define __cmtkQtTriplanarWindow_h_included_

#include <cmtkconfig.h>

#include <qwidget.h>
#include <qslider.h>
#include <qlayout.h>
#include <qpushbutton.h>
#include <qlineedit.h>
#include <qvalidator.h>
#include <qcombobox.h>
#include <qmenubar.h>
#include <qstatusbar.h>
#include <qlabel.h>
#include <qtabwidget.h>
#include <QGridLayout>
#include <QActionGroup>

#include <cmtkMacros.h>
#include <cmtkStudy.h>

#include <cmtkImage.h>
#include <cmtkColormap.h>
#include <cmtkImageToImageRGB.h>

#include <cmtkQtProgress.h>
#include <cmtkQtScrollRenderView.h>
#include <cmtkQtWindowLevelControls.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/** Triplanar image viewer window.
 */
class QtTriplanarWindow :
  /// Inherit from Qt widget.
  public QWidget
{
  Q_OBJECT // we're using slots

public:
  /// Constructor.
  QtTriplanarWindow();

protected:
  /// The Study object we're working on.
  cmtkGetSetMacro(Study::SmartPtr,Study);

  /// Use (linear) interpolation when rendering slices.
  QAction* m_InterpolateAction;

  /// Use (linear) interpolation when rendering slices.
  QAction* m_CrosshairAction;

  /// Fill Null data areas with a checkerbox pattern
  QAction* m_CheckerboxAction;

  /// Zoom factor in percent.
  int m_ZoomFactor;

  /// Update dialog after study change.
  void UpdateDialog();

  /// Batch mode flag: if set, no dialog boxes are shown.
  bool m_BatchMode;

  /// Action group for zoom factor adjustment.
  QActionGroup* m_ZoomActions;

public slots:
  void slotDataChanged( Study::SmartPtr& study );
  void slotColormapChanged( Study::SmartPtr& study );
  void slotSwitchToStudy( Study::SmartPtr& study );
  void slotSwitchToStudyInternal( Study::SmartPtr& study );

protected slots:
  void slotView25();
  void slotView33();
  void slotView50();
  void slotView100();
  void slotView200();
  void slotView300();
  void slotView400();
  void slotView500();

  void slotViewInterpolation();
  void slotViewCrosshair();
  void slotViewCheckerbox();

  void slotExportMenuCmd( int command );
  void slotRenderAll();

  /// Batch mode slots
  void slotSetColormap( const QString& cmap );
  void slotSetWindowLevel( const QString& wl );
  void slotGoToPixel( const QString& xyz );
  void slotGoToLocation( const QString& xyz );
  void slotExportImage( const QString& filename, const int command );
  void slotSetInterpolateMode( const bool mode );
  void slotSetCrosshairMode( const bool mode );
  void slotSetCheckerboardMode( const bool mode );
  void slotSetZoom( const int zoomPercent );

  /// Switch image in axial viewer.
  void slotSwitchImageAx( int imageIndex );
  void slotSwitchImageSa( int imageIndex );
  void slotSwitchImageCo( int imageIndex );

  /// Three-dimensional mouse event.
  void slotMouse3D( Qt::ButtonState, const Vector3D& );

  /// This slot is called when the "Center" button is clicked.
  void slotCenter();

  /// This slot is called when the "Go To Location" button is clicked.
  void slotGoToLocation();

  /// This slot is called when the "Go To Landmark" button is clicked.
  void slotGoToLandmark();

  /// This slot is called when the "Delete Landmark" button is clicked.
  void slotDeleteLandmark();

  /// This slot is called when the "Add Landmark" button is clicked.
  void slotAddLandmark();

  /// This slot is called when the "Export Landmarks" button is clicked.
  void slotExportLandmarks();

  /// This slot is called when the "Import Landmarks" button is clicked.
  void slotImportLandmarks();

protected:
  /// Store volume dimensions here for convenient access.
  unsigned int VolumeDims[3];

  /// The scrolled view we display an image in.
  QtScrollRenderView* ScrollRenderViewAx;
  QtScrollRenderView* ScrollRenderViewSa;
  QtScrollRenderView* ScrollRenderViewCo;

  Image* PipelineImageAx;
  Image* PipelineImageSa;
  Image* PipelineImageCo;

  Colormap* m_Colormap;

  ImageToImageRGB* ImageToImageRGBAx;
  ImageToImageRGB* ImageToImageRGBSa;
  ImageToImageRGB* ImageToImageRGBCo;

  QMenuBar* MenuBar;
  QGridLayout* GridLayout;
  QStatusBar* StatusBar;

  QGridLayout* LandmarksLayout;
  QLineEdit* LocationEntryX;
  QLineEdit* LocationEntryY;
  QLineEdit* LocationEntryZ;
  QDoubleValidator* LocationValidatorX;
  QDoubleValidator* LocationValidatorY;
  QDoubleValidator* LocationValidatorZ;

  QPushButton* GoToLocationButton;
  QPushButton* CenterButton;
  QPushButton* GoToLandmarkButton;
  QPushButton* AddLandmarkButton;
  QPushButton* DeleteLandmarkButton;
  QPushButton* ExportLandmarksButton;
  QPushButton* ImportLandmarksButton;
    
  QComboBox* LandmarkBox;
  QtWindowLevelControls* WindowLevelControls;

  QTabWidget* m_ControlsTab;
  QtProgress* m_ProgressReporter;
 
private:
  /// The pixel grid index of the current location.
  int GridIndex[3];

  /// Status bar output of grid coordinate.
  QLabel* GridIndexInfo;

  /// Update status bar.
  void UpdateGridInfo();
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtTriplanarWindow_h_included_
