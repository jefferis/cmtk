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

#ifndef __cmtkQtScrollRenderView_h_included_
#define __cmtkQtScrollRenderView_h_included_

#include <cmtkconfig.h>

#include <qobject.h>
#include <qscrollarea.h>
#include <qslider.h>
#include <qlabel.h>

#include <QGroupBox>

#include <cmtkQtRenderImageRGB.h>
#include <cmtkImageRGB.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/// Widget that renders an RGB image in a scrolled viewport.
class QtScrollRenderView :
  /// This is essentially a vertical layout of scroll view and index slider.
  public QGroupBox
{
  Q_OBJECT

public:
  /** Constructor.
   * This object will be child of the given parent widget.
   */
  QtScrollRenderView( QWidget *parentWidget, const QString &title = QString::null );

  /// Destructor.
  virtual ~QtScrollRenderView();

  /// Show the optional image index slider at the bottom of the widget.
  void ShowSlider() { this->m_SliderGroupBox->show(); }

  /// Hide the optional image index slider at the bottom of the widget.
  void HideSlider() { this->m_SliderGroupBox->hide(); }

  /// Set left slider label.
  void SetSliderLabelL( const QString& left )
  {
    if ( this->m_LabelL )
      {
      if ( left.isNull() )
	{
	this->m_LabelL->hide();
	}
      else
	{
	this->m_LabelL->setText( left );
	this->m_LabelL->show();
	}
      }
  }

  /// Set right slider label.
  void SetSliderLabelR( const QString& right )
  {
    if ( this->m_LabelR )
      {
      if ( right.isNull() )
	{
	this->m_LabelR->hide();
	}
      else
	{
	this->m_LabelR->setText( right );
	this->m_LabelR->show();
	}
      }
  }

  /// Get render image object.
  QtRenderImageRGB* GetRenderImage() { return RenderImage; }

  /// Get constant render image object.
  const QtRenderImageRGB* GetRenderImage() const { return RenderImage; }

signals:
  /// This signals a change in the index slider's value.
  void indexChanged( int );

  /// This signal is emitted when a mouse button is pressed on the viewer.
  void signalMousePressed( Qt::ButtonState button, int x, int y );

  /// This signal is emitted when a mouse button is pressed on the viewer.
  void signalMouse3D( Qt::ButtonState button, const Vector3D& v );

public slots:
  /// Connect render view to an RGB image object.
  void slotConnectImage( ImageRGB *const image );

  /// Update rendering.
  void slotRender();
  
  /// Set number of slices for the image slider.
  void slotSetNumberOfSlices( unsigned int nSlices );

  /// Set to given slice.
  void slotSetSlice( unsigned int slice );

public:
  /// Get current slice index.
  unsigned int GetSlice() const { return ImageIndexSlider->value(); }

private:
  /// The scrolled viewport.
  QScrollArea* ScrollView;

  /// The actual renderer, if we're rendering to an image object.
  QtRenderImageRGB* RenderImage;

  /// Image index slider.
  QSlider* ImageIndexSlider;

  /// Left slider label.
  QLabel* m_LabelL;

  /// Right slider label.
  QLabel* m_LabelR;

  /// Slider layout.
  QGroupBox* m_SliderGroupBox;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtScrollRenderView_h_included_
