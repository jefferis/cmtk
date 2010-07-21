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

#ifndef __cmtkQtRenderImageRGB_h_included_
#define __cmtkQtRenderImageRGB_h_included_

#include <cmtkconfig.h>

#include "Pipeline/cmtkRenderer.h"

#include <qobject.h>
#include <qwidget.h>
#include <qimage.h>
#include <qpixmap.h>
#include <qpaintdevice.h>

#include "Base/cmtkMacros.h"
#include "Base/cmtkTypes.h"

#include <QMouseEvent>
#include <QPaintEvent>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/** Class to render RGB images in Qt.
 */
class QtRenderImageRGB :
  /// Inherit from QImage for Qt display features.
  public QWidget,
  /// Inherit from Renderer pipeline object.
  public Renderer
{
  Q_OBJECT

  /// Zoom factor.
  igsClassParameter(unsigned int,ZoomFactorPercent);

  /// Flip image in x-direction.
  igsClassParameter(bool,FlipX);

  /// Flip image in y-direction.
  igsClassParameter(bool,FlipY);

  /// Aspect mode constants.
  typedef enum {
    /// No aspect control; render as is.
    AspectNone,
    /// Use x-axis as aspect master.
    AspectX,
    /// Use y-axis as aspect master.
    AspectY,
    /// Automatically select aspect axis.
    AspectAuto
  } AspectMode;

  /// Flip image in y-direction.
  igsClassParameter(AspectMode,ImageAspectMode);

  /// Crosshair mode.
  igsClassParameter(bool,CrosshairMode);

  /// Crosshair position.
  igsClassParameter2Array(unsigned int,CrosshairPosition);

  /// Crosshair colors.
  igsClassParameter2Array(QColor,CrosshairColors);

public:
  /// Constructor.
  QtRenderImageRGB( QWidget *const parent = 0, Qt::WFlags f = 0 );

  /// Destructor.
  virtual ~QtRenderImageRGB();

  /// The actual renderer function.
  virtual void Execute();

  /// Return currently displayed image.
  QPixmap GetPixmap();

  /// Render to given paint device.
  void RenderTo( QPaintDevice *pd );

  /// Render this image.
  virtual void Render()
  {
    this->Renderer::Render();
    this->update();
  }

signals:
  /// This signal is emitted when a mouse button is pressed on the widget.
  void signalMousePressed( Qt::MouseButton button, int x, int y );

  /// This signal is emitted when a mouse button is pressed on the widget.
  void signalMouse3D( Qt::MouseButton button, const Vector3D& v );

protected:
  /// Repaint widget.
  virtual void paintEvent( QPaintEvent *const );

  /// React to mouse clicks (generate a signal).
  virtual void mousePressEvent( QMouseEvent *e );

  /// React to mouse dragging (generate a signal).
  virtual void mouseMoveEvent( QMouseEvent *e );

  /** Capture displayed RGB image.
   * Capture actual display of this widget, potentially with annotations, axes,
   * crosshair, etc.
   */
  virtual ImageRGB* CaptureDisplay();

private:
  /// Intermediate QImage object for painting.
  QImage Image;

  /// Draw crosshair.
  void DrawCrosshair( QPainter &painter, const unsigned int width, const unsigned int height ) const;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtRenderImageRGB_h_included_
