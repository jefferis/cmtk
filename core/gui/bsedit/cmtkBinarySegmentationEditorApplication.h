/*
//
//  Copyright 2010-2011 SRI International
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

#ifndef __cmtkBinarySegmentationEditorApplication_h_included_
#define __cmtkBinarySegmentationEditorApplication_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkXformList.h>

#include <Registration/cmtkReformatVolume.h>

#include <Qt/cmtkQGraphicsPixmapItemEvents.h>

#include <ui_bseditMainWindow.h>

#include <QtGui/QApplication>
#include <QtGui/QMainWindow>
#include <QtGui/QGraphicsView>
#include <QtGui/QGraphicsScene>
#include <QtGui/QGraphicsSceneMouseEvent>
#include <QtGui/QGraphicsLineItem>

namespace
cmtk
{

/// Application class for fusion viewer.
class BinarySegmentationEditorApplication : public QApplication
{
  Q_OBJECT

public:
  /// This class.
  typedef BinarySegmentationEditorApplication Self;

  /// Constructor.
  BinarySegmentationEditorApplication( int& argc, char* argv[] );

private slots:
  /// Update displayed fixed image slice.
  void setFixedSlice( int slice );

  /// Update zoom factor from UI.
  void changeZoom( QAction* action /*!< Action to set new zoom factor. */ );

  /// Change fixed image color map.
  void changeFixedColor( QAction* );

  /// Update slice direction from UI.
  void changeSliceDirection( QAction* action /*!< Action to set new slice direction. */ );

  /// Fixed image black/white has changed.
  void fixedBlackWhiteChanged();

  /// Update slice direction from integer.
  void changeSliceDirection( const int sliceAxis );

  /// Mouse button pressed in one of the graphics views.
  void mousePressed( QGraphicsSceneMouseEvent* event );

private:
  /// Application main window.
  QMainWindow* m_MainWindow;

  /// Designed-generated User Interface for the main window.
  Ui::fviewMainWindow m_MainWindowUI;

  /// Class to bundle data relevant to image objects.
  class Data
  {
  public:
    /// The volume.
    UniformVolume::SmartConstPtr m_Volume;

    /// Data range.
    Types::DataItemRange m_DataRange;

    /// Data for the current image slice.
    UniformVolume::SmartConstPtr m_Slice;

    /// Color map index: this selects one of the predefined color maps.
    int m_ColorMapIndex;

    /// Color table.
    QVector<QRgb> m_ColorTable;

    /// QImage for the current slice.
    QImage m_Image;
  };

  /// The intensity image data.
  Self::Data m_IntensityImage;

  /// Class to bundle data relevant to display objects.
  class Display
  {
  public:
    /// The graphics view (this is a link to the view created from the main window uic).
    QGraphicsView* m_View;
    
    /// The graphics scene for this volume.
    QGraphicsScene* m_Scene;
    
    /// The pixmap graphics item with mouse events.
    QGraphicsPixmapItemEvents* m_PixmapItem;
  };
    
  /// The foreground display data.
  Self::Display m_Foreground;

  /// The foreground display data.
  Self::Display m_Background;

  /// Initialize the view for the given display (fore- or background).
  void InitViewDisplay( Self::Display& display, /*!< Bundled data for given volume.*/ QGraphicsView* view /*!< The view we want to attach this volume to.*/ );

  /// The slice axis (0=x, sagittal; 1=y, coronal; 2=z, axial).
  int m_SliceAxis;

  /// Slice index in the fixed image along the slice axis.
  int m_SliceIndex;

  /// Zoom scale factor.
  float m_ZoomFactor;

  /// Scale factors for non-square pixels.
  FixedVector<3,float> m_ScalePixels;

  /// Update displayed fixed image slice.
  void UpdateImage();

  /// Make a color table based on the color map index.
  void MakeColorTable( Self::Data& data );

  /// Make a QImage from slice data and color table.
  void MakeImage( QImage& image, const UniformVolume& slice, const QVector<QRgb>& colorTable, const float blackLevel, const float whiteLevel );

  /// Update graphics view using a given image.
  void UpdateView( Self::Display& display, QImage& image );

  /// Get 3D coordinate axis corresponding to 2D x axis.
  int GetAxis2DX() const
  {
    static const int idxXtable[3] = { 1, 0, 0 };
    return idxXtable[this->m_SliceAxis];
  }
  
  /// Get 3D coordinate axis corresponding to 2D y axis.
  int GetAxis2DY() const
  {
    static const int idxYtable[3] = { 2, 2, 1 };
    return idxYtable[this->m_SliceAxis];
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkBinarySegmentationEditorApplication_h_included_
