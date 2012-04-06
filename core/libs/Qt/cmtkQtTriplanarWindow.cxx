/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include "cmtkQtTriplanarWindow.h"

#include <Base/cmtkLandmark.h>
#include <Base/cmtkLandmarkList.h>

#include <IO/cmtkClassStream.h>
#include <Qt/cmtkQtIcons.h>

#include <qapplication.h>
#include <qmessagebox.h>
#include <qradiobutton.h>
#include <QMenuBar>
#include <QMenu>
#include <qinputdialog.h>
#include <qfiledialog.h>

#include <QLabel>
#include <QGridLayout>
#include <QPixmap>
#include <QPainter>

#include <iostream>
#include <fstream>
#include <sstream>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtTriplanarWindow::QtTriplanarWindow()
  : QWidget( NULL ),
    m_Study( NULL ),
    m_ZoomFactor( 100 ),
    m_BatchMode( false )
{
  this->setWindowIcon( QtIcons::WindowIcon() );

  this->m_ZoomActions = new QActionGroup( this );
  this->m_ZoomActions->setExclusive( true );
  QAction* action;

  MenuBar = new QMenuBar( this );

  this->ViewMenu = MenuBar->addMenu( "&View" );
  action = ViewMenu->addAction( "25%", this, SLOT( slotView25() ) );
  action->setCheckable( true );
  this->m_ZoomActions->addAction( action );
  action = ViewMenu->addAction( "33%", this, SLOT( slotView33() ) );
  action->setCheckable( true );
  this->m_ZoomActions->addAction( action );
  action = ViewMenu->addAction( "50%", this, SLOT( slotView50() ) );
  action->setCheckable( true );
  this->m_ZoomActions->addAction( action );
  ViewMenu->addSeparator();
  action = ViewMenu->addAction( "100%", this, SLOT( slotView100() ) );
  action->setCheckable( true );
  action->setChecked( true );  
  this->m_ZoomActions->addAction( action );
  ViewMenu->addSeparator();
  action = ViewMenu->addAction( "200%", this, SLOT( slotView200() ) );
  action->setCheckable( true );
  this->m_ZoomActions->addAction( action );
  action = ViewMenu->addAction( "300%", this, SLOT( slotView300() ) );
  action->setCheckable( true );
  this->m_ZoomActions->addAction( action );
  action = ViewMenu->addAction( "400%", this, SLOT( slotView400() ) );
  action->setCheckable( true );
  this->m_ZoomActions->addAction( action );
  action = ViewMenu->addAction( "500%", this, SLOT( slotView500() ) );
  action->setCheckable( true );
  this->m_ZoomActions->addAction( action );
  ViewMenu->addSeparator();
  (this->m_InterpolateAction = ViewMenu->addAction( "&Interpolation", this, SLOT( slotViewInterpolation() ) ))->setCheckable( true );
  this->m_InterpolateAction->setChecked( false );
  (this->m_CrosshairAction = ViewMenu->addAction( "&Crosshair", this, SLOT( slotViewCrosshair() ) ))->setCheckable( true );
  this->m_CrosshairAction->setChecked( true );
  (this->m_CheckerboxAction = ViewMenu->addAction( "C&heckerbox", this, SLOT( slotViewCheckerbox() ) ))->setCheckable( true );
  this->m_CheckerboxAction->setChecked( true );
  
  this->ExportMenu = MenuBar->addMenu( "&Export" );
  ExportMenu->addAction( "&Axial" )->setData( QVariant( 1 ) );
  ExportMenu->addAction( "&Coronal" )->setData( QVariant( 2 ) );
  ExportMenu->addAction( "&Sagittal" )->setData( QVariant( 3 ) );
  ExportMenu->addSeparator();
  ExportMenu->addAction( "&Panel")->setData( QVariant( 4 ) );
  QObject::connect( ExportMenu, SIGNAL( triggered( QAction* ) ), this, SLOT( slotExportMenuCmd( QAction* ) ) );

  MenuBar->show();

  StatusBar = new QStatusBar( this );
  StatusBar->show();
  GridIndexInfo = new QLabel( StatusBar );
  GridIndex[0] = GridIndex[1] = GridIndex[2] = 0;
  StatusBar->addWidget( GridIndexInfo, 1 );

  GridLayout = new QGridLayout( this );
  GridLayout->setMenuBar( MenuBar );
  GridLayout->addWidget( StatusBar, 2, 0, 1, 2 );

  ScrollRenderViewAx = new QtScrollRenderView( this, "Axial" );
  ScrollRenderViewAx->SetSliderLabelL( "I" );
  ScrollRenderViewAx->SetSliderLabelR( "S" );
  ScrollRenderViewAx->GetRenderImage()->SetFlipX( true );
  ScrollRenderViewAx->GetRenderImage()->SetFlipY( true );
  ScrollRenderViewAx->ShowSlider();
  ScrollRenderViewAx->GetRenderImage()->SetCrosshairMode( this->m_CrosshairAction->isChecked() );
  ScrollRenderViewAx->GetRenderImage()->SetCrosshairColors( QColor( 255,0,0 ), QColor( 0,255,0 ) );
  GridLayout->addWidget( ScrollRenderViewAx, 1, 0 );

  ScrollRenderViewSa = new QtScrollRenderView( this, "Sagittal" );
  ScrollRenderViewSa->ShowSlider();
  ScrollRenderViewSa->SetSliderLabelL( "L" );
  ScrollRenderViewSa->SetSliderLabelR( "R" );
  ScrollRenderViewSa->GetRenderImage()->SetCrosshairMode( this->m_CrosshairAction->isChecked() );
  ScrollRenderViewSa->GetRenderImage()->SetCrosshairColors( QColor( 0,255,0 ), QColor( 0,0,255 ) );
  GridLayout->addWidget( ScrollRenderViewSa, 0, 1 );

  ScrollRenderViewCo = new QtScrollRenderView( this, "Coronal" );
  ScrollRenderViewCo->ShowSlider();
  ScrollRenderViewCo->SetSliderLabelL( "P" );
  ScrollRenderViewCo->SetSliderLabelR( "A" );
  ScrollRenderViewCo->GetRenderImage()->SetFlipX( true );
  ScrollRenderViewCo->GetRenderImage()->SetCrosshairMode( this->m_CrosshairAction->isChecked() );
  ScrollRenderViewCo->GetRenderImage()->SetCrosshairColors( QColor( 255,0,0 ), QColor( 0,0,255 ) );
  GridLayout->addWidget( ScrollRenderViewCo, 0, 0 );

  this->m_ControlsTab = new QTabWidget( this );
  this->m_ControlsTab->setContentsMargins( 10, 10, 10, 10 );
  GridLayout->addWidget( this->m_ControlsTab, 1, 1 );

  QWidget *landmarksTab = new QWidget( this->m_ControlsTab );
  this->m_ControlsTab->addTab( landmarksTab, "Landmarks" );

  LandmarksLayout = new QGridLayout( landmarksTab );

  LocationEntryX = new QLineEdit( landmarksTab );
  LocationValidatorX = new QDoubleValidator( LocationEntryX );
  LocationValidatorX->setDecimals( 6 );
  LocationEntryX->setValidator( LocationValidatorX );
  LandmarksLayout->addWidget( LocationEntryX, 0, 0 );

  QObject::connect( LocationEntryX, SIGNAL( returnPressed() ), this, SLOT( slotGoToLocation() ) );

  LocationEntryY = new QLineEdit( landmarksTab );
  LocationValidatorY = new QDoubleValidator( LocationEntryY );
  LocationEntryY->setValidator( LocationValidatorY );
  LocationValidatorY->setDecimals( 6 );
  LandmarksLayout->addWidget( LocationEntryY, 0, 1 );

  QObject::connect( LocationEntryY, SIGNAL( returnPressed() ), this, SLOT( slotGoToLocation() ) );

  LocationEntryZ = new QLineEdit( landmarksTab );
  LocationValidatorZ = new QDoubleValidator( LocationEntryZ );
  LocationEntryZ->setValidator( LocationValidatorZ );
  LocationValidatorZ->setDecimals( 6 );
  LandmarksLayout->addWidget( LocationEntryZ, 0, 2 );

  QObject::connect( LocationEntryZ, SIGNAL( returnPressed() ), this, SLOT( slotGoToLocation() ) );

  CenterButton = new QPushButton( landmarksTab );
  CenterButton->setText( "Center" );
  LandmarksLayout->addWidget( CenterButton, 1, 1 );
  QObject::connect( CenterButton, SIGNAL( clicked() ), this, SLOT( slotCenter() ) );

  GoToLocationButton = new QPushButton( landmarksTab );
  GoToLocationButton->setText( "Go To" );
  LandmarksLayout->addWidget( GoToLocationButton, 1, 2 );
  QObject::connect( GoToLocationButton, SIGNAL( clicked() ), this, SLOT( slotGoToLocation() ) );

  LandmarkBox = new QComboBox( landmarksTab );
  LandmarksLayout->addWidget( LandmarkBox, 2, 0, 1, 2 );
  LandmarkBox->setEnabled( false );
  QObject::connect( LandmarkBox, SIGNAL( currentIndexChanged(int) ), this, SLOT( slotGoToLandmark() ) );

  GoToLandmarkButton = new QPushButton( landmarksTab );
  GoToLandmarkButton->setText( "Go To" );
  LandmarksLayout->addWidget( GoToLandmarkButton, 3, 0 );
  QObject::connect( GoToLandmarkButton, SIGNAL( clicked() ), this, SLOT( slotGoToLandmark() ) );
  GoToLandmarkButton->setEnabled( false );

  DeleteLandmarkButton = new QPushButton( landmarksTab );
  DeleteLandmarkButton->setText( "Delete" );
  LandmarksLayout->addWidget( DeleteLandmarkButton, 3, 1 );
  QObject::connect( DeleteLandmarkButton, SIGNAL( clicked() ), this, SLOT( slotDeleteLandmark() ) );
  DeleteLandmarkButton->setEnabled( false );
  
  AddLandmarkButton = new QPushButton( landmarksTab );
  AddLandmarkButton->setText( "Add..." );
  LandmarksLayout->addWidget( AddLandmarkButton, 3, 2 );
  QObject::connect( AddLandmarkButton, SIGNAL( clicked() ), this, SLOT( slotAddLandmark() ) );
  
  ExportLandmarksButton = new QPushButton( landmarksTab );
  ExportLandmarksButton->setText( "Export Landmarks..." );
  LandmarksLayout->addWidget( ExportLandmarksButton, 5, 0 );
  QObject::connect( ExportLandmarksButton, SIGNAL( clicked() ), this, SLOT( slotExportLandmarks() ) );
  ExportLandmarksButton->setEnabled( false );

  ImportLandmarksButton = new QPushButton( landmarksTab );
  ImportLandmarksButton->setText( "Import Landmarks..." );
  LandmarksLayout->addWidget( ImportLandmarksButton, 5, 2 );
  QObject::connect( ImportLandmarksButton, SIGNAL( clicked() ), this, SLOT( slotImportLandmarks() ) );
  ImportLandmarksButton->setEnabled( true );

  // create the window/level tab
  WindowLevelControls = new QtWindowLevelControls( this->m_ControlsTab );
  QObject::connect( WindowLevelControls, SIGNAL( colormap( Study::SmartPtr& ) ), this, SLOT( slotColormapChanged( Study::SmartPtr& ) ) );
  this->m_ControlsTab->addTab( WindowLevelControls, "Window/Level" );

  // create the visualization pipeline
  this->m_Colormap = Colormap::New();
  this->m_Colormap->SetStandardColormap( PALETTE_GRAY );
  this->m_Colormap->SetDataRange( 0, 2000 );

  // set up axial viewer
  PipelineImageAx = Image::New();
  ImageToImageRGBAx = ImageToImageRGB::New();
  ImageToImageRGBAx->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  ImageToImageRGBAx->SetInput( PipelineImageAx );
  ImageToImageRGBAx->SetColormap( this->m_Colormap );
  ScrollRenderViewAx->slotConnectImage( ImageToImageRGBAx->GetOutput() );

  QObject::connect( ScrollRenderViewAx, SIGNAL( indexChanged( int ) ), this, SLOT( slotSwitchImageAx( int ) ) );
  QObject::connect( ScrollRenderViewAx, SIGNAL( signalMouse3D( Qt::MouseButton, const Vector3D& ) ), this, SLOT( slotMouseAx( Qt::MouseButton, const Vector3D& ) ) );
  
  // set up sagittal viewer
  PipelineImageSa = Image::New();
  ImageToImageRGBSa = ImageToImageRGB::New();
  ImageToImageRGBSa->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  ImageToImageRGBSa->SetInput( PipelineImageSa );
  ImageToImageRGBSa->SetColormap( this->m_Colormap );
  ScrollRenderViewSa->slotConnectImage( ImageToImageRGBSa->GetOutput() );

  PipelineImageCo = Image::New();
  QObject::connect( ScrollRenderViewSa, SIGNAL( indexChanged( int ) ), this, SLOT( slotSwitchImageSa( int ) ) );
  QObject::connect( ScrollRenderViewSa, SIGNAL( signalMouse3D( Qt::MouseButton, const Vector3D& ) ), this, SLOT( slotMouseSa( Qt::MouseButton, const Vector3D& ) ) );

  // set up coronal viewer
  ImageToImageRGBCo = ImageToImageRGB::New();
  ImageToImageRGBCo->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  ImageToImageRGBCo->SetInput( PipelineImageCo );
  ImageToImageRGBCo->SetColormap( this->m_Colormap );
  ScrollRenderViewCo->slotConnectImage( ImageToImageRGBCo->GetOutput() );

  QObject::connect( ScrollRenderViewCo, SIGNAL( indexChanged( int ) ), this, SLOT( slotSwitchImageCo( int ) ) );
  QObject::connect( ScrollRenderViewCo, SIGNAL( signalMouse3D( Qt::MouseButton, const Vector3D& ) ), this, SLOT( slotMouseCo( Qt::MouseButton, const Vector3D& ) ) );

  this->m_ProgressReporter = new QtProgress( this );
}


void
QtTriplanarWindow::slotSetZoom( const int zoomPercent )
{
  this->m_ZoomFactor = zoomPercent;
  ScrollRenderViewAx->GetRenderImage()->SetZoomFactorPercent( zoomPercent );
  ScrollRenderViewCo->GetRenderImage()->SetZoomFactorPercent( zoomPercent );
  ScrollRenderViewSa->GetRenderImage()->SetZoomFactorPercent( zoomPercent );
  this->slotRenderAll();
}

void
QtTriplanarWindow::slotSetCheckerboardMode( const bool mode )
{
  this->m_CheckerboxAction->setChecked( mode );
  this->slotViewCheckerbox();
}
      
void
QtTriplanarWindow::slotSetCrosshairMode( const bool mode )
{
  this->m_CrosshairAction->setChecked( mode );
  this->slotViewCrosshair();
}
      
void
QtTriplanarWindow::slotSetInterpolateMode( const bool mode )
{
  this->m_InterpolateAction->setChecked( mode );
  this->slotViewInterpolation();
}
      
void
QtTriplanarWindow::slotExportMenuCmd( QAction* action )
{
  const int mode = action->data().toInt();

  QString title( "Choose filename" );
  switch ( mode ) 
    {
    case 1: // Axial image.
      title = "Axial image export";
      break;
    case 2: // Coronal image.
      title = "Coronal image export";
      break;
    case 3: // Sagittal.
      title = "Sagittal image export";
      break;
    case 4: 
      title = "Panel image export";
      break;
    }
  
  QString filename( "image.png" );
  filename = QFileDialog::getSaveFileName( this, title, filename, "Portable Network Graphic (*.png);; Tagged Image File Format (*.tif);; Portable Pixmap (*.ppm *.pgm);; JPEG (*.jpg)" );
  
  if ( !filename.isEmpty() ) 
    {
    this->slotExportImage( filename, mode );
    }
}


void
QtTriplanarWindow::slotExportImage( const QString& filename, const int command )
{  
  QPixmap pixmap;
  switch ( command ) 
    {
    case 1: // Axial image.
      pixmap = ScrollRenderViewAx->GetRenderImage()->GetPixmap();
      break;
    case 2: // Coronal image.
      pixmap = ScrollRenderViewCo->GetRenderImage()->GetPixmap();
      break;
    case 3: // Sagittal.
      pixmap = ScrollRenderViewSa->GetRenderImage()->GetPixmap();
      break;
    case 4: 
    { // Panel.
    QPixmap pixmapAx = ScrollRenderViewAx->GetRenderImage()->GetPixmap();
    QPixmap pixmapSa = ScrollRenderViewSa->GetRenderImage()->GetPixmap();
    QPixmap pixmapCo = ScrollRenderViewCo->GetRenderImage()->GetPixmap();
    
    pixmap = QPixmap( pixmapCo.width() + pixmapSa.width(), pixmapCo.height() + pixmapAx.height() );
    QPainter painter( &pixmap );

    painter.drawPixmap( 0, 0, pixmapCo.width(), pixmapCo.height(), pixmapCo );
    painter.drawPixmap( pixmapCo.width(), 0, pixmapSa.width(), pixmapSa.height(), pixmapSa );
    painter.drawPixmap( 0, pixmapCo.height(), pixmapAx.width(), pixmapAx.height(), pixmapAx );
    painter.fillRect( pixmapCo.width(), pixmapCo.height(), pixmapSa.width(), pixmapAx.height(), Qt::black );
    break;
    }
    }
  
  QString format = filename.section( ".", -1 ).toUpper();
  if ( format.isEmpty() )
    format = "PNG";
  if ( !pixmap.save( filename, format.toLatin1() ) )
    {
    if ( this->m_BatchMode )
      std::cerr << "WARNING: saving file failed." << std::endl;
    else
      QMessageBox::warning( this, "Save failed", "Error saving file" );
    }
}

void
QtTriplanarWindow::slotRenderAll()
{
  ScrollRenderViewAx->slotRender();
  ScrollRenderViewCo->slotRender();
  ScrollRenderViewSa->slotRender();
}

void 
QtTriplanarWindow::slotSwitchToStudy( Study::SmartPtr& study )
{
  this->m_Study = study;
  if ( this->m_Study ) 
    {
    qApp->setOverrideCursor( Qt::WaitCursor );
    this->m_Study->ReadVolume(  true /* reload */, AnatomicalOrientation::ORIENTATION_STANDARD );
    qApp->restoreOverrideCursor();
    
    if ( !this->m_BatchMode )
      {
      while ( ! this->m_Study->GetVolume() ) 
	{
	int button = QMessageBox::warning ( NULL, "Error", "Could not read image data for this study.", QMessageBox::Retry, QMessageBox::Abort );
	if ( button == QMessageBox::Abort ) break;
	}
      }

    if ( this->m_Study->GetVolume() ) 
      {
      this->SetStudy( this->m_Study );
      this->WindowLevelControls->slotSetStudy( this->m_Study );
      this->slotCenter();
      this->slotColormapChanged( this->m_Study );
      this->UpdateDialog();
      this->show();
      }
    else
      {
      if ( this->m_BatchMode )
	StdErr << "ERROR: could not read image " << this->m_Study->GetFileSystemPath() << "\n";
      }
    
    // if study has landmarks, put them into combo box.
    LandmarkBox->clear();
    const LandmarkList* ll = this->m_Study->GetLandmarkList();
    if ( ll ) 
      {
      LandmarkList::const_iterator ll_it = ll->begin();
      while ( ll_it != ll->end() ) 
	{
	LandmarkBox->addItem( ll_it->m_Name.c_str() );
	++ll_it;
	}
      }
    LandmarkBox->setEnabled( LandmarkBox->count() );
    GoToLandmarkButton->setEnabled( LandmarkBox->count() );
    DeleteLandmarkButton->setEnabled( LandmarkBox->count() );
    ExportLandmarksButton->setEnabled( LandmarkBox->count() );
  }
}

void 
QtTriplanarWindow::slotSwitchToStudyInternal( Study::SmartPtr& study )
{
  this->m_Study = study;
  if ( this->m_Study ) 
    {
    this->m_Study->ReadVolume( false /* reload */, AnatomicalOrientation::ORIENTATION_STANDARD );
    
    while ( ! this->m_Study->GetVolume() ) 
      {
      int button = QMessageBox::warning ( NULL, "Error", "Could not read image data for this study.", QMessageBox::Retry, QMessageBox::Abort );
      if ( button == QMessageBox::Abort ) break;
      }

    if ( this->m_Study->GetVolume() ) 
      {
      this->SetStudy( this->m_Study );
      this->WindowLevelControls->slotSetStudy( this->m_Study );

      this->slotSwitchImageAx( ScrollRenderViewAx->GetSlice() );
      this->slotSwitchImageSa( ScrollRenderViewSa->GetSlice() );
      this->slotSwitchImageCo( ScrollRenderViewCo->GetSlice() );
      
      this->UpdateDialog();
      this->show();
      }
    }
}

void
QtTriplanarWindow::UpdateDialog()
{
  if ( this->m_Study ) 
    {
    const UniformVolume *volume = this->m_Study->GetVolume();
    if ( volume ) 
      {
      VolumeDims = volume->GetDims();
      
      ScrollRenderViewAx->slotSetNumberOfSlices( VolumeDims[AXIS_Z] );
      ScrollRenderViewSa->slotSetNumberOfSlices( VolumeDims[AXIS_X] );
      ScrollRenderViewCo->slotSetNumberOfSlices( VolumeDims[AXIS_Y] );
      
      LocationValidatorX->setBottom( 0 );
      LocationValidatorX->setTop( volume->Size[AXIS_X] );
      LocationValidatorY->setBottom( 0 );
      LocationValidatorY->setTop( volume->Size[AXIS_Y] );
      LocationValidatorZ->setBottom( 0 );
      LocationValidatorZ->setTop( volume->Size[AXIS_Z] );
      } 
    else
      {
      qWarning( "QtTriplanarWindow::UpdateDialog called with no image data loaded.\n" );
      }
    if ( this->m_CrosshairAction->isChecked() ) 
      {
      this->slotRenderAll();
      }
    
    QString caption;
    this->setWindowTitle( caption.sprintf( "CMTK Triplanar Viewer: %s", this->m_Study->GetName() ) );
    this->show();
    }
}

void
QtTriplanarWindow::slotDataChanged( Study::SmartPtr& study )
{
  if ( study != this->m_Study ) return;
  this->slotGoToLocation();
}

void
QtTriplanarWindow::slotSwitchImageAx( int imageIndex )
{
  const UniformVolume *volume = this->m_Study->GetVolume();

  if ( volume ) 
    {
    ScalarImage::SmartPtr sliceImage( volume->GetOrthoSlice( AXIS_Z, imageIndex ) );
    
    if ( sliceImage ) 
      {
      if ( ! this->m_CheckerboxAction->isChecked() )
	sliceImage->GetPixelData()->ReplacePaddingData( 0.0 );
      
      sliceImage->AdjustToIsotropic( volume->GetMinDelta(), this->m_InterpolateAction->isChecked() );
      PipelineImageAx->SetFromScalarImage( sliceImage );
      }
    sliceImage = ScalarImage::SmartPtr::Null();

    LocationEntryZ->setText( QString::number( volume->GetPlaneCoord( AXIS_Z, imageIndex ) ) );
    GridIndex[2] = imageIndex;
    this->UpdateGridInfo();
    
    if ( this->m_CrosshairAction->isChecked() ) 
      {
      this->slotGoToLocation();
      } 
    else
      {
      ScrollRenderViewAx->slotRender();
      }
    } 
  else
    {
    qWarning( "QtTriplanarWindow::slotSwitchImageAx called with no image data loaded.\n" );
    }
}

void
QtTriplanarWindow::slotSwitchImageSa( int imageIndex )
{
  const UniformVolume *volume = this->m_Study->GetVolume();

  if ( volume ) 
    {
    ScalarImage *sliceImage = volume->GetOrthoSlice( AXIS_X, imageIndex );
    
    if ( sliceImage ) 
      {
      if ( ! this->m_CheckerboxAction->isChecked() )
	sliceImage->GetPixelData()->ReplacePaddingData( 0.0 );
      
      sliceImage->Mirror( false /* horizontal */, true /* vertical */ );
      sliceImage->AdjustToIsotropic( volume->GetMinDelta(), this->m_InterpolateAction->isChecked() );
      PipelineImageSa->SetFromScalarImage( sliceImage );
      delete sliceImage;
      }
    
    LocationEntryX->setText( QString::number( volume->GetPlaneCoord( AXIS_X, imageIndex ) ) );
    GridIndex[0] = imageIndex;
    this->UpdateGridInfo();
    
    if ( this->m_CrosshairAction->isChecked() ) 
      {
      this->slotGoToLocation();
      } 
    else
      {
      ScrollRenderViewSa->slotRender();
      }
    } 
  else
    {
    qWarning( "QtTriplanarWindow::slotSwitchImageSa called with no image data loaded.\n" );
    }
}

void
QtTriplanarWindow::slotSwitchImageCo( int imageIndex )
{
  const UniformVolume *volume = this->m_Study->GetVolume();

  if ( volume ) 
    {
    ScalarImage *sliceImage = volume->GetOrthoSlice( AXIS_Y, imageIndex );
    
    if ( sliceImage ) 
      {
      if ( ! this->m_CheckerboxAction->isChecked() )
	sliceImage->GetPixelData()->ReplacePaddingData( 0.0 );
      
      sliceImage->Mirror( false /* horizontal */, true /* vertical */ );
      sliceImage->AdjustToIsotropic( volume->GetMinDelta(), this->m_InterpolateAction->isChecked() );
      PipelineImageCo->SetFromScalarImage( sliceImage );
      delete sliceImage;
      }
    
    LocationEntryY->setText( QString::number( volume->GetPlaneCoord( AXIS_Y, imageIndex ) ) );
    GridIndex[1] = imageIndex;
    this->UpdateGridInfo();
    
    if ( this->m_CrosshairAction->isChecked() ) 
      {
      this->slotGoToLocation();
      } 
    else
      {
      ScrollRenderViewCo->slotRender();
      }
    
    } 
  else
    {
    qWarning( "QtTriplanarWindow::slotSwitchImageCo called with no image data loaded.\n" );
    }
}

void
QtTriplanarWindow
::slotGoToPixel( const QString& xyz )
{
  int x, y, z;
  if ( 3 != sscanf( xyz.toLatin1(), "%d,%d,%d", &x, &y, &z ) )
    {
    qWarning( "QtTriplanarWindow::slotGoToPixel needs pixel index as 'x,y,z'.\n" );
    }
  else
    {
    this->slotSwitchImageSa( x );
    this->slotSwitchImageCo( y );
    this->slotSwitchImageAx( z );
    }
}

void
QtTriplanarWindow::slotGoToLocation( const QString& xyz )
{
  float v[3];
  if ( 3 != sscanf( xyz.toLatin1(), "%f,%f,%f", v, v+1, v+2 ) )
    {
    qWarning( "QtTriplanarWindow::slotGoToLocation needs 3D coordinate as 'x,y,z'.\n" );
    }
  else
    {
    this->slotMouse3D( Qt::LeftButton, UniformVolume::CoordinateVectorType( v ) );
    }
}

void
QtTriplanarWindow
::slotSetColormap( const QString& cmap )
{
  for ( unsigned int colormapIndex = 0; Colormap::StandardColormaps[colormapIndex]; ++colormapIndex ) 
    {
    if ( cmap == QString( Colormap::StandardColormaps[colormapIndex] ) )
      {
      this->m_Colormap->SetStandardColormap( colormapIndex );
      this->slotRenderAll();
      break;
      }
    }
}

void
QtTriplanarWindow
::slotSetWindowLevel( const QString& wl )
{
  float window, level;
  if ( 2 != sscanf( wl.toLatin1(), "%f:%f", &window, &level ) )
    {
    qWarning( "QtTriplanarWindow::slotSetWindowLevel needs 'window:level'.\n" );
    }
  else
    {
    this->m_Colormap->SetDataRange( level-0.5*window, level+0.5*window );
    this->slotRenderAll();
    }
}

void 
QtTriplanarWindow
::slotMouse3D( Qt::MouseButton, const Vector3D& v )
{
// if we don't have a study yet, simply ignore.
  if ( ! this->m_Study )
    return;

  const UniformVolume *volume = this->m_Study->GetVolume();

  unsigned int i=0, j=0;
  PipelineImageAx->ProjectPixel( v, i, j );
  ScrollRenderViewAx->GetRenderImage()->SetCrosshairPosition( i, j );

  PipelineImageSa->ProjectPixel( v, i, j ); 
  ScrollRenderViewSa->GetRenderImage()->SetCrosshairPosition( i, j );

  PipelineImageCo->ProjectPixel( v, i, j ); 
  ScrollRenderViewCo->GetRenderImage()->SetCrosshairPosition( i, j );

  if ( volume )
    {
    const unsigned int sliceSa = volume->GetClosestCoordIndex( AXIS_X, v[AXIS_X] );
    ScrollRenderViewSa->slotSetSlice( sliceSa );
    ScrollRenderViewSa->slotRender();
    
    const unsigned int sliceCo = volume->GetClosestCoordIndex( AXIS_Y, v[AXIS_Y] );
    ScrollRenderViewCo->slotSetSlice( sliceCo );
    ScrollRenderViewCo->slotRender();
    
    const unsigned int sliceAx = volume->GetClosestCoordIndex( AXIS_Z, v[AXIS_Z] );
    ScrollRenderViewAx->slotSetSlice( sliceAx );
    ScrollRenderViewAx->slotRender();
    }
}

void 
QtTriplanarWindow
::slotMouseAx( Qt::MouseButton, const Vector3D& v )
{
// if we don't have a study yet, simply ignore.
  if ( ! this->m_Study )
    return;

  const UniformVolume *volume = this->m_Study->GetVolume();

  unsigned int i=0, j=0;
  PipelineImageAx->ProjectPixel( v, i, j );
  ScrollRenderViewAx->GetRenderImage()->SetCrosshairPosition( i, j );

  if ( volume )
    {
    const unsigned int sliceSa = volume->GetClosestCoordIndex( AXIS_X, v[AXIS_X] );
    ScrollRenderViewSa->slotSetSlice( sliceSa );
    ScrollRenderViewSa->slotRender();
    
    const unsigned int sliceCo = volume->GetClosestCoordIndex( AXIS_Y, v[AXIS_Y] );
    ScrollRenderViewCo->slotSetSlice( sliceCo );
    ScrollRenderViewCo->slotRender();
    }
}

void 
QtTriplanarWindow
::slotMouseSa( Qt::MouseButton, const Vector3D& v )
{
// if we don't have a study yet, simply ignore.
  if ( ! this->m_Study )
    return;

  const UniformVolume *volume = this->m_Study->GetVolume();

  unsigned int i=0, j=0;
  PipelineImageSa->ProjectPixel( v, i, j );
  ScrollRenderViewSa->GetRenderImage()->SetCrosshairPosition( i, j );

  if ( volume )
    {
    const unsigned int sliceAx = volume->GetClosestCoordIndex( AXIS_Z, v[AXIS_Z] );
    ScrollRenderViewAx->slotSetSlice( sliceAx );
    ScrollRenderViewAx->slotRender();
    
    const unsigned int sliceCo = volume->GetClosestCoordIndex( AXIS_Y, v[AXIS_Y] );
    ScrollRenderViewCo->slotSetSlice( sliceCo );
    ScrollRenderViewCo->slotRender();
    }
}

void 
QtTriplanarWindow
::slotMouseCo( Qt::MouseButton, const Vector3D& v )
{
// if we don't have a study yet, simply ignore.
  if ( ! this->m_Study )
    return;

  const UniformVolume *volume = this->m_Study->GetVolume();

  unsigned int i=0, j=0;
  PipelineImageCo->ProjectPixel( v, i, j );
  ScrollRenderViewCo->GetRenderImage()->SetCrosshairPosition( i, j );

  if ( volume )
    {
    const unsigned int sliceAx = volume->GetClosestCoordIndex( AXIS_Z, v[AXIS_Z] );
    ScrollRenderViewAx->slotSetSlice( sliceAx );
    ScrollRenderViewAx->slotRender();
    
    const unsigned int sliceSa = volume->GetClosestCoordIndex( AXIS_X, v[AXIS_X] );
    ScrollRenderViewSa->slotSetSlice( sliceSa );
    ScrollRenderViewSa->slotRender();
    }
}

void
QtTriplanarWindow::slotCenter()
{
  const UniformVolume *volume = this->m_Study->GetVolume();
  if ( ! volume ) return;

  // Pretend there was a button event at the given location
  this->slotMouse3D( Qt::LeftButton, volume->GetCenterCropRegion() );
}

void
QtTriplanarWindow::slotGoToLocation()
{
  const UniformVolume *volume = this->m_Study->GetVolume();
  if ( ! volume ) return;

  // Pretend there was a button event at the given location
  const double location[3] = { LocationEntryX->text().toDouble(), LocationEntryY->text().toDouble(), LocationEntryZ->text().toDouble() };
  this->slotMouse3D( Qt::LeftButton, UniformVolume::CoordinateVectorType( location ) );
}

void
QtTriplanarWindow::slotGoToLandmark()
{
  if ( ! this->m_Study ) return;

  const LandmarkList *ll = this->m_Study->GetLandmarkList();
  if ( ! ll ) return;

  LandmarkList::ConstIterator lm = ll->FindByName( LandmarkBox->currentText().toStdString() );
  if ( lm != ll->end() ) 
    {
    this->slotMouse3D( Qt::LeftButton, lm->m_Location );
    }
}

void
QtTriplanarWindow::slotDeleteLandmark()
{
}

void
QtTriplanarWindow::slotExportLandmarks()
{
  if ( ! this->m_Study ) return;

  LandmarkList::SmartPtr ll = this->m_Study->GetLandmarkList();
  if ( ! ll ) return;

  QString path = QFileDialog::getSaveFileName( this, "Save Landmarks File" );
  
  if ( ! path.isEmpty() ) 
    {
    std::ofstream stream( path.toLatin1().constData() );
    
    if ( stream.good() )
      {
      for ( LandmarkList::ConstIterator it = ll->begin(); it != ll->end(); ++it )
	{
	stream << it->m_Location[0] << "\t" << it->m_Location[1] << "\t" << it->m_Location[2] << "\t" << it->m_Name << std::endl;
	}
      stream.close();
      } 
    else 
      {
      QMessageBox::critical( NULL, "Error", "Could not open file for writing.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
      }
    }
}

void
QtTriplanarWindow::slotImportLandmarks()
{
  if ( ! this->m_Study ) return;

  LandmarkList::SmartPtr ll = this->m_Study->GetLandmarkList();
  if ( ! ll ) 
    {
    ll = LandmarkList::SmartPtr( new LandmarkList );
    this->m_Study->SetLandmarkList( ll );
    }
  
  QString path = QFileDialog::getOpenFileName( this, "Open Landmarks File", QString(), "All Files (*.*)" );
  
  if ( ! path.isEmpty() ) 
    {
    std::ifstream stream( path.toLatin1().constData() );

    unsigned int cnt = 0;
    if ( stream.good() ) 
      {
      while ( ! stream.eof() )
	{
	Landmark::SpaceVectorType xyz;
	stream >> xyz[0] >> xyz[1] >> xyz[2];

	char name[128];
	stream.getline( name, 128, '\n' );

	if ( ! strlen(name) )
	  {
	  sprintf( name, "LM-%04d", cnt++ );
	  }

	ll->push_back( Landmark( name, xyz ) );
	LandmarkBox->addItem( name );
	}
      
      LandmarkList::ConstIterator lm = ll->begin();
      if ( lm != ll->end() )
	{
	this->LandmarkBox->setCurrentIndex( this->LandmarkBox->findText( lm->m_Name.c_str() ) );
	this->slotMouse3D( Qt::LeftButton, lm->m_Location );
	}
      
      LandmarkBox->setEnabled( true );
      GoToLandmarkButton->setEnabled( true );
      DeleteLandmarkButton->setEnabled( true );
      ExportLandmarksButton->setEnabled( true );
      
      stream.close();
      } 
    else 
      {
      QMessageBox::critical( NULL, "Error", "Could not open file for reading.", QMessageBox::Ok, Qt::NoButton, Qt::NoButton );
      }
    }
}

void
QtTriplanarWindow::slotAddLandmark()
{
  if ( ! this->m_Study ) return;

  LandmarkList::SmartPtr ll = this->m_Study->GetLandmarkList();
  if ( ! ll ) 
    {
    ll = LandmarkList::SmartPtr( new LandmarkList );
    this->m_Study->SetLandmarkList( ll );
    }
  
  bool ok;
  QString name = QInputDialog::getText( this, "Add New Landmark", "Enter new landmark name:", QLineEdit::Normal, QString::null, &ok );
  if ( ok && !name.isEmpty() ) 
    {
    Types::Coordinate location[3] = { LocationEntryX->text().toDouble(), LocationEntryY->text().toDouble(), LocationEntryZ->text().toDouble() };
    ll->push_back( Landmark( name.toStdString(), Landmark::SpaceVectorType( location ) ) );
    LandmarkBox->addItem( name );
    LandmarkBox->setCurrentIndex( this->LandmarkBox->findText( name ) );
    
    LandmarkBox->setEnabled( true );
    GoToLandmarkButton->setEnabled( true );
    DeleteLandmarkButton->setEnabled( true );
    ExportLandmarksButton->setEnabled( true );
    }
}

void
QtTriplanarWindow::slotColormapChanged( Study::SmartPtr& study )
{
  if ( this->m_Study && (this->m_Study == study) ) 
    {
    this->m_Colormap->SetFromStudy( this->m_Study.GetPtr() );
    this->slotRenderAll();
    }
}

void
QtTriplanarWindow::UpdateGridInfo()
{
  if ( ! (this->m_Study && this->m_Study->GetVolume()) ) return;

  QString str = "OUTSIDE";

  const UniformVolume* volume = this->m_Study->GetVolume();
  if ( volume->IndexIsInRange( GridIndex[0], GridIndex[1], GridIndex[2] ) )
    {
    const FixedVector<3,float> v = volume->IndexToPhysical( GridIndex );

    Types::DataItem value;
    if ( volume->GetDataAt( value, GridIndex[0], GridIndex[1], GridIndex[2] ) )
      str.sprintf( "Pixel Index: [%d,%d,%d] RAS: [%g,%g,%g] Value: %g", GridIndex[0], GridIndex[1], GridIndex[2], v[0], v[1], v[2], value );
    else
      str.sprintf( "Pixel Index: [%d,%d,%d] RAS: [%g,%g,%g]", GridIndex[0], GridIndex[1], GridIndex[2], v[0], v[1], v[2] );
    }
  GridIndexInfo->setText( str );
}

} // namespace cmtk
