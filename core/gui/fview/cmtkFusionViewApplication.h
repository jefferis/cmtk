/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkFusionViewApplication_h_included_
#define __cmtkFusionViewApplication_h_included_

#include <cmtkconfig.h>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkXformList.h>

#include <Registration/cmtkReformatVolume.h>

#include <Qt/cmtkQGraphicsPixmapItemEvents.h>

#include <ui_fviewMainWindow.h>

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
class FusionViewApplication : public QApplication
{
  Q_OBJECT

public:
  /// This class.
  typedef FusionViewApplication Self;

  /// Constructor.
  FusionViewApplication( int argc, char* argv[] );

private slots:
  /// Update displayed fixed image slice.
  void setFixedSlice( int slice );

  /// Update moving image transparency.
  void setTransparency( int slice );

  /// Update flag for displaying linked cursors.
  void setLinkedCursorFlag( bool flag );

  /// Update flag for affine-only transformations.
  void setAffineOnly( bool affineOnly );

  /// Update zoom factor from UI.
  void changeZoom( QAction* action /*!< Action to set new zoom factor. */ );

  /// Update interpolator from UI.
  void changeInterpolator( QAction* action /*!< Action to set new interpolator. */ );

  /// Update slice direction from UI.
  void changeSliceDirection( QAction* action /*!< Action to set new slice direction. */ );

  /// Fixed image black/white has changed.
  void fixedBlackWhiteChanged();

  /// Moving image black/white has changed.
  void movingBlackWhiteChanged();

  /// Update slice direction from integer.
  void changeSliceDirection( const int sliceAxis );

  /// Mouse button pressed in one of the graphics views.
  void mousePressed( QGraphicsSceneMouseEvent* event );

private:
  /// Application main window.
  QMainWindow* m_MainWindow;

  /// Designed-generated User Interface for the main window.
  Ui::fviewMainWindow m_MainWindowUI;

  /// Class to bundle fixed and moving image objects.
  class Data
  {
  public:
    /// The volume.
    UniformVolume::SmartConstPtr m_Volume;

    /// Data range.
    Types::DataItemRange m_DataRange;

    /// Data for the current image slice.
    UniformVolume::SmartConstPtr m_Slice;

    /// Color table.
    QVector<QRgb> m_ColorTable;

    /// QImage for the current slice.
    QImage m_Image;

    /// The graphics view (this is a link to the view created from the main window uic).
    QGraphicsView* m_View;

    /// The graphics scene for this volume.
    QGraphicsScene* m_Scene;

    /// The pixmap graphics item with mouse events.
    QGraphicsPixmapItemEvents* m_PixmapItem;

    /// The line items for the cross cursor.
    QGraphicsLineItem* m_CursorLines[2];
  };

  /// The fixed volume data.
  Self::Data m_Fixed;

  /// The fixed volume data.
  Self::Data m_Moving;

  /// Initialize the view data for the given volume (fixed or moving).
  void InitViewData( Self::Data& data, /*!< Bundled data for given volume.*/ QGraphicsView* view /*!< The view we want to attach this volume to.*/ );

  /// The list of concatenated transformations.
  XformList m_XformList;

  /** The list of all-affine concatenated transformations.
   * these are the same transformations as in m_XformList, but every nonrigid
   * transformation therein is replaced with its affine initializer.
   */
  XformList m_XformListAllAffine;

  /// Flag to apply only affine transformation components.
  bool m_AffineOnly;

  /// The slice axis (0=x, sagittal; 1=y, coronal; 2=z, axial).
  int m_SliceAxis;

  /// Slice index in the fixed image along the slice axis.
  int m_SliceIndex;

  /// Interpolator for the moving image.
  Interpolators::InterpolationEnum m_Interpolator;

  /// QImage for the current fused slice.
  QImage m_FusedImage;

  /// Zoom scale factor.
  float m_ZoomFactor;

  /// Scale factors for non-square pixels.
  FixedVector<3,float> m_ScalePixels;

  /// Moving image transparency.
  float m_Transparency;

  /// Flag for linked cursor display.
  float m_CursorDisplayed;

  /// Linked cursor position in 3D.
  FixedVector<3,float> m_CursorPosition;

  /// Update displayed fixed image slice.
  void UpdateFixedImage();

  /// Update interpolated moving image slice.
  void UpdateMovingSlice();

  /// Update displayed moving image slice.
  void UpdateMovingImage();

  /// Make a QImage from slice data and color table.
  void MakeImage( QImage& image, const UniformVolume& slice, const QVector<QRgb>& colorTable, const float blackLevel, const float whiteLevel );

  /// Update graphics view using a given image.
  void UpdateView( Self::Data& data, QImage& image );

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

#endif // #ifndef __cmtkFusionViewApplication_h_included_
