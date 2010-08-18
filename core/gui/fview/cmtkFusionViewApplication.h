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

#include <QtGui/QApplication>
#include <QtGui/QMainWindow>

#include <Base/cmtkUniformVolume.h>
#include <Base/cmtkXformList.h>

#include <Registration/cmtkReformatVolume.h>

#include <ui_fviewMainWindow.h>

namespace
cmtk
{

/// Application class for fusion viewer.
class FusionViewApplication : public QApplication
{
  Q_OBJECT

public:
  /// Constructor.
  FusionViewApplication( int argc, char* argv[] );

private slots:
  /// Update displayed fixed image slice.
  void setFixedSlice( int slice );

  /// Update moving image transparency.
  void setTransparency( int slice );

  /// Update moving image transparency.
  void setLinkedCursorFlag( bool flag );

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

private:
  /// Application main window.
  QMainWindow* m_MainWindow;

  /// Designed-generated User Interface for the main window.
  Ui::fviewMainWindow m_MainWindowUI;

  /// The fixed volume.
  UniformVolume::SmartConstPtr m_FixedVolume;

  /// The moving volume.
  UniformVolume::SmartConstPtr m_MovingVolume;

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

  /// Data for the current fixed image slice.
  UniformVolume::SmartConstPtr m_FixedSlice;

  /// Data for the current moving image slice.
  UniformVolume::SmartConstPtr m_MovingSlice;

  /// Color table for fixed image.
  QVector<QRgb> m_ColorTableFix;

  /// Color table for moving image.
  QVector<QRgb> m_ColorTableMov;

  /// QImage for the current fixed slice.
  QImage m_FixedImage;

  /// QImage for the current moving slice.
  QImage m_MovingImage;

  /// QImage for the current fused slice.
  QImage m_FusedImage;

  /// Zoom scale factor.
  float m_ZoomFactor;

  /// Scale factors for non-square pixels.
  float m_ScalePixels[2];

  /// Moving image transparency.
  float m_Transparency;

  /// Flag for linked cursor display.
  float m_CursorDisplayed;

  /// Linked cursor position.
  float m_CursorPosition[2];

  /// Update displayed fixed image slice.
  void UpdateFixedImage();

  /// Update interpolated moving image slice.
  void UpdateMovingSlice();

  /// Update displayed moving image slice.
  void UpdateMovingImage();

  /// Make a QImage from slice data and color table.
  void MakeImage( QImage& image, const UniformVolume& slice, const QVector<QRgb>& colorTable, const float blackLevel, const float whiteLevel );

  /// Update graphics view using slice data, black and white levels.
  void UpdateView( QGraphicsView* view, const QImage& image );
};

} // namespace cmtk

#endif // #ifndef __cmtkFusionViewApplication_h_included_
