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

#include <cmtkQtStudyWidget.h>

#include <qapplication.h>
#include <qmessagebox.h>
#include <qcombobox.h>
#include <qlayout.h>
#include <QHBoxLayout>
#include <QVBoxLayout>

namespace
cmtk
{

QtStudyWidget::QtStudyWidget
( QWidget* parent, Qt::WFlags flags )
  : QWidget( parent, flags ),
    m_Study( NULL )
{
  // create the visualization pipeline
  PipelineImage = Image::New();
  this->m_Colormap = Colormap::New();
  this->m_Colormap->SetStandardColormap( PALETTE_GRAY );
  this->m_Colormap->SetDataRange( 0, 2000 );
  this->m_ImageToImageRGB = ImageToImageRGB::New();
  this->m_ImageToImageRGB->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  this->m_ImageToImageRGB->SetInput( PipelineImage );
  this->m_ImageToImageRGB->SetColormap( this->m_Colormap );

  ScrollRenderView = new QtScrollRenderView( this );
  ScrollRenderView->ShowSlider();
  ScrollRenderView->slotConnectImage( this->m_ImageToImageRGB->GetOutput() );

  QObject::connect( ScrollRenderView, SIGNAL( indexChanged( int ) ), this, SLOT( slotSwitchImage( int ) ) );

  QBoxLayout* topLevelLayout = new QHBoxLayout( this );
  topLevelLayout->addWidget( ScrollRenderView );

  QBoxLayout* rightLayout = new QVBoxLayout;
  topLevelLayout->addLayout( rightLayout, 0 );

  WindowLevelBox = new QtWindowLevelControls( this );
  QObject::connect( WindowLevelBox, SIGNAL( colormap( Study::SmartPtr& ) ), SIGNAL( colormap( Study::SmartPtr& ) ) );
  rightLayout->addWidget( WindowLevelBox );

  QObject::connect( this, SIGNAL( colormap( Study::SmartPtr& ) ), SLOT( slotUpdateColormap() ) );
}

void 
QtStudyWidget::slotSetStudy( Study::SmartPtr& study )
{
  this->m_Study = study;

  if ( this->m_Study ) 
    {
    qApp->setOverrideCursor( Qt::WaitCursor );
    this->m_Study->ReadVolume( false /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
    qApp->restoreOverrideCursor();
    
    Volume::SmartPtr volume( NULL );
    while ( ! (volume = this->m_Study->GetVolume()) ) 
      {
      int button = QMessageBox::warning ( NULL, "Error", "Could not read image data for this study.", QMessageBox::Retry, QMessageBox::Abort );
      if ( button == QMessageBox::Abort ) break;
      this->m_Study->ReadVolume( false /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
      }
    
    if ( volume ) 
      {
      WindowLevelBox->slotSetStudy( this->m_Study );
      ScrollRenderView->slotSetNumberOfSlices( volume->GetDims( AXIS_Z ) );
      ScrollRenderView->slotRender();
      this->show();
      }
    }
}

void
QtStudyWidget::slotDataChanged( Study::SmartPtr& study )
{
  if ( this->m_Study == study ) 
    {
    this->WindowLevelBox->slotSetStudy( this->m_Study );
    this->slotSwitchImage( this->ImageIndex );
    this->ScrollRenderView->slotRender();
    }
}

void 
QtStudyWidget::slotSwitchImage( int imageIndex )
{
  if ( ! this->m_Study ) return;
  
  Volume::SmartPtr volume = this->m_Study->GetVolume();
  if ( volume ) 
    {
    this->ImageIndex = imageIndex;
    
    ScalarImage::SmartPtr sliceImage( volume->GetOrthoSlice( AXIS_Z, this->ImageIndex ) );
    if ( sliceImage )
      {
      sliceImage->AdjustAspect();
      PipelineImage->SetFromScalarImage( sliceImage.GetPtr() );
      }
    ScrollRenderView->slotRender();
    } 
  else 
    {
    qWarning( "QtStudyWidget::slotSwitchImage called with no image data loaded.\n" );
    }
}

void
QtStudyWidget::slotUpdateColormap()
{
  if ( this->m_Study )
    {
    this->m_Colormap->SetFromStudy( this->m_Study );
    ScrollRenderView->slotRender();
    }
}

} // namespace cmtk
