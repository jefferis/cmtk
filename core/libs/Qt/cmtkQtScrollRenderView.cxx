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

#include "cmtkQtScrollRenderView.h"

#include <qlabel.h>
#include <QFrame>
#include <QGridLayout>
#include <QVBoxLayout>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtScrollRenderView::QtScrollRenderView
( QWidget *parentWidget, const QString& title )
  : QGroupBox( parentWidget ),
    RenderImage( NULL )
{
  if ( ! parentWidget )
    qFatal( "No parent widget in QtScrollRenderView constructor." );

  if ( title != QString::null ) 
    {
    this->setAlignment ( Qt::AlignLeft );
    this->setTitle( title );
    } 
  
  ScrollView = new QScrollArea( this ); //, "ScrollView", Qt::WResizeNoErase|Qt::WStaticContents  );
  RenderImage = new QtRenderImageRGB( this );
  ScrollView->setWidget( RenderImage );
  ScrollView->setFrameStyle( QFrame::NoFrame );

  // export viewer's mouse signal
  QObject::connect( RenderImage,SIGNAL(signalMousePressed( Qt::MouseButton, int, int )), SIGNAL(signalMousePressed( Qt::MouseButton, int, int )) );
  QObject::connect( RenderImage, SIGNAL(signalMouse3D( Qt::MouseButton, const Vector3D& )), SIGNAL(signalMouse3D( Qt::MouseButton, const Vector3D& )) );
  
  // lock minimum size so we can display 256x256 with not scrollbars.
  RenderImage->setMinimumSize( 256, 256 );

  this->m_SliderGroupBox = new QGroupBox( this );
  this->m_SliderGroupBox->hide();

  QGridLayout* sliderBoxLayout = new QGridLayout( this->m_SliderGroupBox );
  sliderBoxLayout->setContentsMargins( 0,0,0,0 );
  
  this->ImageIndexSlider = new QSlider( this->m_SliderGroupBox  );
  this->ImageIndexSlider->setOrientation( Qt::Horizontal );
  this->ImageIndexSlider->setDisabled( true );
  sliderBoxLayout->addWidget( this->ImageIndexSlider, 0, 1 );

  this->m_LabelL = new QLabel( this->m_SliderGroupBox );
  sliderBoxLayout->addWidget( this->m_LabelL, 0, 0 );
  this->m_LabelR = new QLabel( this->m_SliderGroupBox );
  sliderBoxLayout->addWidget( this->m_LabelR, 0, 2 );

  QVBoxLayout *vbox = new QVBoxLayout;  
  vbox->setContentsMargins( 0,0,0,0 );
  vbox->addWidget( this->ScrollView );
  vbox->addWidget( this->m_SliderGroupBox );
  vbox->setSpacing( 0 );
  this->setLayout(vbox);

  // export slider's signal
//  QObject::connect( ImageIndexSlider, SIGNAL( valueChanged( int ) ), SIGNAL( indexChanged( int ) ) );
}

QtScrollRenderView::~QtScrollRenderView()
{
}

void QtScrollRenderView::slotConnectImage( ImageRGB *const image )
{
  if ( RenderImage ) 
    {
    RenderImage->SetInput( image );
    } 
  else 
    {
    qWarning( "RenderImage is NULL in QtScrollRenderView::ConnectRenderView." );
    }
}

void QtScrollRenderView::slotRender()
{
  if ( RenderImage ) 
    {
    RenderImage->Render();
    } 
  else
    {
    qWarning( "RenderImage is NULL in QtScrollRenderView::Render." );
    }
}

void QtScrollRenderView::slotSetNumberOfSlices( unsigned int nSlices )
{
  if ( nSlices ) 
    {
    ImageIndexSlider->setEnabled( true );
    ImageIndexSlider->setMinimum( 0 );
    ImageIndexSlider->setMaximum( nSlices - 1 );

    if ( ImageIndexSlider->value() < 0 || ImageIndexSlider->value() >= (int)nSlices )
      {
      ImageIndexSlider->setValue( nSlices / 2 );
      }
    this->ImageIndexSlider->setDisabled( false );
    } 
  else
    {
    ImageIndexSlider->setDisabled( true );
    }
}

void QtScrollRenderView::slotSetSlice( unsigned int slice )
{
  if ( slice <= static_cast<unsigned>(ImageIndexSlider->maximum()) ) 
    {
    ImageIndexSlider->setValue( slice );
    }
}

} // namespace cmtk
