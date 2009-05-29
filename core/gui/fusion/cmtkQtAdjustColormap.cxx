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

#include <cmtkQtAdjustColormap.h>

#include <cmtkQtFusionGlobal.h>

#include <qapplication.h>
#include <qmessagebox.h>
#include <qlayout.h>
#include <q3scrollview.h>
#include <q3vgroupbox.h>
#include <qradiobutton.h>
#include <qmenubar.h>

#include <QEvent>
#include <Q3HBoxLayout>
#include <Q3VBoxLayout>
#include <Q3PopupMenu>
#include <Q3ButtonGroup>

namespace
cmtk
{

QtAdjustColormap::QtAdjustColormap()
  : QWidget( NULL, "Adjust Colormap" ),
    m_Study( NULL )
{
  this->setIcon( QtFusionGlobal::WindowIcon() );

  Q3PopupMenu* ViewMenu = new Q3PopupMenu();
  ViewMenu->insertItem( "&100%" );
  ViewMenu->insertItem( "&200%" );
  ViewMenu->insertItem( "&300%" );
  ViewMenu->insertItem( "&400%" );
  
  QMenuBar* MenuBar = new QMenuBar( this );
  MenuBar->insertItem( "&View", ViewMenu );
  MenuBar->show();

  Q3BoxLayout* topLevelLayout = new Q3HBoxLayout( this );
  topLevelLayout->setMenuBar( MenuBar );

  ScrollRenderView = new QtScrollRenderView( this );
  ScrollRenderView->ShowSlider();
  topLevelLayout->addWidget( ScrollRenderView );

  Q3BoxLayout* rightLayout = new Q3VBoxLayout( topLevelLayout );
  WindowLevelBox = new QtWindowLevelControls( this );
  rightLayout->addWidget( WindowLevelBox );
  QObject::connect( WindowLevelBox, SIGNAL( signalChanged(float,float,float) ), this, SLOT( slotWindowChanged() ) );

  // create palette selection buttons.
  Q3ButtonGroup *paletteGroup = new Q3HButtonGroup( this );
  paletteGroup->setTitle( "Colormap" );
  const char** colorMapNames = Colormap::StandardColormaps;
  int paletteIdx = 0;
  while ( colorMapNames[paletteIdx] ) 
    {
    paletteGroup->insert( new QRadioButton( colorMapNames[paletteIdx], paletteGroup ), paletteIdx );
    ++paletteIdx;
    }
  paletteGroup->setRadioButtonExclusive( true );
  QObject::connect( paletteGroup, SIGNAL( pressed( int ) ), this, SLOT( slotSetStandardColormap( int ) ) );
  rightLayout->addWidget( paletteGroup );

  // insert stretch below all widgets in the right column; left has 
  // auto-stretching scroll view.
  rightLayout->insertStretch( -1 );
  rightLayout->activate();
  rightLayout->setResizeMode( QLayout::SetFixedSize );

  // create the visualization pipeline
  PipelineImage = Image::New();
  this->m_Colormap = Colormap::New();
  this->m_Colormap->SetStandardColormap( PALETTE_GRAY );
  this->m_Colormap->SetDataRange( 0, 2000 );
  this->m_ImageToImageRGB = ImageToImageRGB::New();
  this->m_ImageToImageRGB->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  this->m_ImageToImageRGB->SetInput( PipelineImage );
  this->m_ImageToImageRGB->SetColormap( this->m_Colormap );
  ScrollRenderView->slotConnectImage( this->m_ImageToImageRGB->GetOutput() );

  QObject::connect( ScrollRenderView, SIGNAL( indexChanged( int ) ), this, SLOT( slotSwitchImage( int ) ) );
}

QtAdjustColormap::~QtAdjustColormap()
{
  qWarning( "QtAdjustColormap destroyed.\n" );
}

void 
QtAdjustColormap::SwitchToStudy( Study::SmartPtr& study )
{
  if ( this->m_Study )
    emit newIcon( this->m_Study, ScrollRenderView->GetRenderImage()->GetIcon( 24 ) );

  if ( study ) 
    {
    qApp->setOverrideCursor( Qt::waitCursor );
    study->ReadVolume( false /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
    qApp->restoreOverrideCursor();
    
    while ( ! study->GetVolume() ) 
      {
      int button = QMessageBox::warning ( NULL, "Error", "Could not read image data for this study.", QMessageBox::Retry, QMessageBox::Abort );
      if ( button == QMessageBox::Abort ) break;
      }
    if ( study->GetVolume() ) 
      {
      this->SetStudy( study );
      WindowLevelBox->slotSetStudy( this->m_Study );
      this->UpdateDialog();
      this->show();
      emit newIcon( this->m_Study, ScrollRenderView->GetRenderImage()->GetIcon( 24 ) );
      }
    }
}

void
QtAdjustColormap::UpdateDialog()
{
  if ( this->m_Study ) 
    {
    const Volume *volume = this->m_Study->GetVolume();
    if ( volume ) 
      {
      ScrollRenderView->slotSetNumberOfSlices( volume->GetDims( AXIS_Z ) );
      } 
    else
      {
      qWarning( "QtAdjustColormap::UpdateDialog called with no image data loaded.\n" );
      }
    ScrollRenderView->slotRender();
    
    QString caption;
    this->setCaption( caption.sprintf( "Adjust Colormap: %s", this->m_Study->GetName() ) );
    this->show();
    }
}

void
QtAdjustColormap::slotSwitchImage( int imageIndex )
{
  const Volume *volume = this->m_Study->GetVolume();
  
  if ( volume ) 
    {
    ScalarImage::SmartPtr sliceImage( volume->GetOrthoSlice( AXIS_Z, imageIndex ) );
    if ( sliceImage ) 
      {
      sliceImage->AdjustAspect();
      PipelineImage->SetFromScalarImage( sliceImage );
      }
    ScrollRenderView->slotRender();
    } 
  else 
    {
    qWarning( "QtAdjustColormap::SwitchImage called with no image data loaded.\n" );
    }
}

void
QtAdjustColormap::slotWindowChanged()
{
  if ( this->m_Study ) 
    {
    this->m_Colormap->SetDataRange( this->m_Study->GetBlack(), this->m_Study->GetWhite() );
    //    this->m_Colormap->SetReverse( this->m_Study->GetReverseColormap() );
    this->m_Colormap->SetGamma( this->m_Study->GetGamma() );
    ScrollRenderView->slotRender();
    emit colormapChanged( this->m_Study );
    }
}

void 
QtAdjustColormap::slotSetStandardColormap( int paletteIdx )
{
  if ( this->m_Study ) 
    {
    this->m_Study->SetStandardColormap( paletteIdx );
    this->m_Colormap->SetStandardColormap( paletteIdx );
    ScrollRenderView->slotRender();
    emit colormapChanged( this->m_Study );
    }
}

void
QtAdjustColormap::leaveEvent( QEvent *event )
{
  this->QWidget::leaveEvent( event );
  // if we lost the focus, let's update the study's icon in the main window.
  emit newIcon( this->m_Study, ScrollRenderView->GetRenderImage()->GetIcon( 24 ) );
}

} // namespace cmtk
