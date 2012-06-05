/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkQtRenderImageRGB.h"

#include <qpainter.h>
#include <qpixmap.h>
#include <qcursor.h>

#include <QPaintEvent>
#include <QMouseEvent>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtRenderImageRGB::QtRenderImageRGB
( QWidget *const parent, Qt::WFlags f )
  : QWidget( parent, f ),
    ZoomFactorPercent( 100 ),
    FlipX( false ), FlipY( false ),
    ImageAspectMode( AspectNone ),
    CrosshairMode( false ),
    Image()
{
  CrosshairPosition[0] = CrosshairPosition[1] = 0;

  this->setBaseSize( 512, 512 );
  //  this->setMouseTracking( true );
  this->setCursor( QCursor( Qt::CrossCursor ) );
}

QtRenderImageRGB::~QtRenderImageRGB()
{
}

void
QtRenderImageRGB::paintEvent( QPaintEvent *const )
{
  // update modified time to force rendering
  this->UpdateModifiedTime();
  this->Update();
  this->RenderTo( this );
  this->UpdateExecuteTime();
}

void
QtRenderImageRGB::Execute()
{
  if ( Input == NULL ) return;
}

void
QtRenderImageRGB::RenderTo( QPaintDevice *pd )
{
  if ( ! Input ) return;

  if ( Input->GetAlphaChannel() != IMAGE_RGBA ) 
    {
    return;
    }
  
  unsigned char* imageDataRGB = (unsigned char*) Input->GetDataPtr();
  if ( imageDataRGB == NULL ) return;

  unsigned int width, height;
  Input->GetDims( width, height );

  this->setFixedSize( ZoomFactorPercent * width / 100, ZoomFactorPercent * height / 100 );
  Image = QImage( width, height, QImage::Format_RGB32 );
  memcpy( Image.bits(), imageDataRGB, width * height * 4 );

  if ( FlipX || FlipY )
    Image = Image.mirrored( FlipX, FlipY );

  if ( ZoomFactorPercent != 100 ) 
    {
    Image = Image.scaled( width * ZoomFactorPercent / 100, height * ZoomFactorPercent / 100 );
    }
  
  QPainter painter( pd );
  painter.drawImage( 0, 0, Image );
  
  if ( CrosshairMode )
    this->DrawCrosshair( painter, width, height );
}

void
QtRenderImageRGB::DrawCrosshair
( QPainter &painter, const unsigned int width, const unsigned int height ) 
  const
{
  int crosshairX = ( FlipX ) ? width-1-CrosshairPosition[0] : CrosshairPosition[0];
  crosshairX = static_cast<int>( (0.5+crosshairX) * ZoomFactorPercent / 100 );

  int crosshairY = ( FlipY ) ? height-1-CrosshairPosition[1] : CrosshairPosition[1];
  crosshairY = static_cast<int>( (0.5+crosshairY) * ZoomFactorPercent / 100 );

  const int realWidth = static_cast<int>( 1.0 * width * ZoomFactorPercent / 100 );
  const int realHeight = static_cast<int>( 1.0 * height * ZoomFactorPercent / 100 );
  
  painter.setPen( CrosshairColors[0] );
  painter.drawLine( 0, crosshairY, realWidth-1, crosshairY );
  
  painter.setPen( CrosshairColors[1] );
  painter.drawLine( crosshairX, 0, crosshairX, realHeight-1 );
}

QPixmap
QtRenderImageRGB::GetPixmap()
{
  if ( Input == NULL ) return QPixmap();
  
  QPixmap capture( ZoomFactorPercent * Input->GetDims( AXIS_X ) / 100, ZoomFactorPercent * Input->GetDims( AXIS_Y ) / 100 );
  this->RenderTo( &capture );
  return capture;
}

/// React to mouse clicks (generate a signal).
void 
QtRenderImageRGB
::mousePressEvent( QMouseEvent *e ) 
{
  unsigned int scaledX = (e->x()-static_cast<int>(ZoomFactorPercent/200)) * 100 / ZoomFactorPercent;
  if ( Input && this->FlipX )
    {
    scaledX = Input->GetDims( AXIS_X ) - 1 - scaledX;
    }

  unsigned int scaledY = (e->y()-static_cast<int>(ZoomFactorPercent/200)) * 100 / ZoomFactorPercent;
  if ( Input && this->FlipY )
    {
    scaledY = Input->GetDims( AXIS_Y ) - 1 - scaledY;
    }
  
  emit signalMousePressed( e->button(), scaledX, scaledY );
  Vector3D v;
  Input->GetPixelLocation( v, scaledX, scaledY );
  emit signalMouse3D( e->button(), v );
  e->accept();
}

void
QtRenderImageRGB
::mouseMoveEvent( QMouseEvent *e ) 
{
  unsigned int scaledX = (e->x()-static_cast<int>(ZoomFactorPercent/200)) * 100 / ZoomFactorPercent;
  if ( Input && this->FlipX )
    {
    scaledX = Input->GetDims( AXIS_X ) - 1 - scaledX;
    }

  unsigned int scaledY = (e->y()-static_cast<int>(ZoomFactorPercent/200)) * 100 / ZoomFactorPercent;
  if ( Input && this->FlipY )
    {
    scaledY = Input->GetDims( AXIS_Y ) - 1 - scaledY;
    }  
  
  emit signalMousePressed( e->button(), scaledX, scaledY );
  Vector3D v;
  Input->GetPixelLocation( v, scaledX, scaledY );

  emit signalMouse3D( e->button(), v );
  e->accept();
}

} // namespace cmtk
