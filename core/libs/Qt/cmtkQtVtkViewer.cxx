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

#include <cmtkQtVtkViewer.h>

#if defined(CMTK_HAVE_QT) && defined(CMTK_HAVE_VTK)

#include <qpushbutton.h>

#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>

#include <cmtkQtIcons.h>
#include <cmtkMathUtil.h>

#include <Q3HBoxLayout>
#include <Q3GridLayout>
#include <Q3PopupMenu>
#include <Q3VBoxLayout>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtVtkViewer::QtVtkViewer
( QWidget *const parent, const char* name, Qt::WFlags flags ) 
  : QWidget( parent, name, flags ),
    Interactor( NULL ),
    WindowToImage( NULL ),
    TIFFWriter( NULL )
{
  Renderer[0][0] = Renderer[0][1] = Renderer[1][0] = Renderer[1][1];

  RenderWindow = NULL;
  this->Construct();
}

QtVtkViewer::QtVtkViewer
( const int nrows, const int ncols, 
  QWidget *const parent, const char* name, Qt::WFlags flags )
  : QWidget( parent, name, flags ),
    Interactor( NULL ),
    WindowToImage( NULL ),
    TIFFWriter( NULL )
{
  Renderer[0][0] = Renderer[0][1] = Renderer[1][0] = Renderer[1][1];

  RenderWindow = NULL;
  this->Construct( nrows, ncols );
}

void
QtVtkViewer::Construct( const int nrows, const int ncols )
{
  this->setCaption( "VTK Viewer" );
  this->setIcon( QtIcons::WindowIcon() );
  
  MenuBar = new QMenuBar( this );

  Q3PopupMenu* viewMenu = new Q3PopupMenu;
  Q3PopupMenu* viewLayoutMenu = new Q3PopupMenu;
  viewLayoutMenu->insertItem( "&Single", this, SLOT( slotSetLayoutSingle() ) );
  viewLayoutMenu->insertItem( "&Panel", this, SLOT( slotSetLayoutPanel() ) );
  viewMenu->insertItem( "&Layout", viewLayoutMenu );
  viewMenu->insertSeparator();
  viewMenu->insertItem( "&Close", this, SLOT( close() ) );

  MenuBar->insertItem( "View", viewMenu );

  //  MenuBar->insertItem( "Scene", sceneMenu );
  MenuBar->show();
  
  MasterLayout = new Q3GridLayout( this, nrows, ncols );
  MasterLayout->setMenuBar( MenuBar );
  
  Q3BoxLayout* ButtonLayoutL = new Q3VBoxLayout();
  MasterLayout->addLayout( ButtonLayoutL, 0, 0 );

  QPushButton* button = new QPushButton( this );
  button->setText( "0" );
  ButtonLayoutL->insertStretch( 2, 1 );
  button->setFixedSize( 32, 32 );
  ButtonLayoutL->addWidget( button );    
  QObject::connect( button, SIGNAL( clicked() ), this, SLOT( slotResetCamera() ) );

  ButtonLayoutL->insertSpacing( 1, 32 );

  const char* directionTxt[] = { "FH", "HF", "LR", "RL", "AP", "PA" };
  for ( unsigned int direction = 0; direction < 6; ++direction ) 
    {
    QPushButton* button = new QPushButton( this );
    button->setText( directionTxt[direction] );
    button->setFixedSize( 32, 32 );
    ButtonLayoutL->addWidget( button );    
    switch ( direction ) 
      {
      case 0:
      default:
	QObject::connect( button, SIGNAL( clicked() ), this, SLOT( slotSetCameraDirectionFH() ) );
	break;
      case 1:
	QObject::connect( button, SIGNAL( clicked() ), this, SLOT( slotSetCameraDirectionHF() ) );
	break;
      case 2:
	QObject::connect( button, SIGNAL( clicked() ), this, SLOT( slotSetCameraDirectionLR() ) );
	break;
      case 3:
	QObject::connect( button, SIGNAL( clicked() ), this, SLOT( slotSetCameraDirectionRL() ) );
	break;
      case 4:
	QObject::connect( button, SIGNAL( clicked() ), this, SLOT( slotSetCameraDirectionAP() ) );
	break;
      case 5:
	QObject::connect( button, SIGNAL( clicked() ), this, SLOT( slotSetCameraDirectionPA() ) );
	break;
      }
    }

  QPushButton* dollyPlus10 = new QPushButton( this );
  dollyPlus10->setText( "+10" );
  dollyPlus10->setFixedSize( 32, 32 );
  QObject::connect( dollyPlus10, SIGNAL( clicked() ), this, SLOT( slotDollyPlus10() ) );
  ButtonLayoutL->addWidget( dollyPlus10 );

  QPushButton* dollyMinus10 = new QPushButton( this );
  dollyMinus10->setText( "-10" );
  dollyMinus10->setFixedSize( 32, 32 );
  QObject::connect( dollyMinus10, SIGNAL( clicked() ), this, SLOT( slotDollyMinus10() ) );
  ButtonLayoutL->addWidget( dollyMinus10 );

  ButtonLayoutL->insertStretch( 8, 1 );

  Q3BoxLayout* ButtonLayoutR = new Q3VBoxLayout();
  MasterLayout->addLayout( ButtonLayoutR, 0, 2 );

  QPushButton* elevationPlus10 = new QPushButton( this );
  elevationPlus10->setText( "+10" );
  elevationPlus10->setFixedSize( 32, 32 );
  QObject::connect( elevationPlus10, SIGNAL( clicked() ), this, SLOT( slotElevationPlus10() ) );
  ButtonLayoutR->addWidget( elevationPlus10 );

  QPushButton* elevationMinus10 = new QPushButton( this );
  elevationMinus10->setText( "-10" );
  elevationMinus10->setFixedSize( 32, 32 );
  QObject::connect( elevationMinus10, SIGNAL( clicked() ), this, SLOT( slotElevationMinus10() ) );
  ButtonLayoutR->addWidget( elevationMinus10 );

  ButtonLayoutR->insertStretch( 0, 1 );

  Q3BoxLayout* ButtonLayoutB = new Q3HBoxLayout();
  MasterLayout->addLayout( ButtonLayoutB, 1, 1 );

  QPushButton* azimuthPlus10 = new QPushButton( this );
  azimuthPlus10->setText( "+10" );
  azimuthPlus10->setFixedSize( 32, 32 );
  QObject::connect( azimuthPlus10, SIGNAL( clicked() ), this, SLOT( slotAzimuthPlus10() ) );
  ButtonLayoutB->addWidget( azimuthPlus10 );

  QPushButton* azimuthMinus10 = new QPushButton( this );
  azimuthMinus10->setText( "-10" );
  azimuthMinus10->setFixedSize( 32, 32 );
  QObject::connect( azimuthMinus10, SIGNAL( clicked() ), this, SLOT( slotAzimuthMinus10() ) );
  ButtonLayoutB->addWidget( azimuthMinus10 );

  ButtonLayoutB->insertStretch( 0, 1 );

  RenderWindow = new vtkQtRenderWindow( this );
  for ( int row = 0; row < 2; ++row ) 
    {
    for ( int col = 0; col < 2; ++col ) 
      {
      Renderer[row][col] = vtkRenderer::New();
      Renderer[row][col]->SetBackground(0.1, 0.2, 0.4);
      if ( row && col )
	Renderer[row][col]->InteractiveOn();
      else
	Renderer[row][col]->InteractiveOff();
      RenderWindow->AddRenderer( Renderer[row][col] );
      }
    }
  MasterLayout->addWidget( RenderWindow, 0, 1 );

  const double vn00[3] = { 0, 0, -1 };
  const double vup00[3] = { 0, 1, 0 };
  this->SetCamera( vn00, vup00, 0, 0 );

  const double vn01[3] = { 1, 0, 0 };
  const double vup01[3] = { 0, 0, 1 };
  this->SetCamera( vn01, vup01, 0, 1 );

  const double vn10[3] = { 0, 1, 0 };
  const double vup10[3] = { 0, 0, 1 };
  this->SetCamera( vn10, vup10, 1, 0 );

  Interactor = vtkQtRenderWindowInteractor::New();
  Interactor->SetRenderWindow( RenderWindow );
}

void 
QtVtkViewer::slotResetCamera()
{  
  for ( int row = 0; row < 2; ++row )
    for ( int col = 0; col < 2; ++col ) 
      {
      vtkCamera* camera = Renderer[row][col]->GetActiveCamera();
      
      const double* vn = camera->GetViewPlaneNormal(); 
      const double* vup = camera->GetViewUp();
      this->SetCamera( vn, vup, row, col );
      }
  
  this->slotUpdateRenderers();
}

void 
QtVtkViewer::slotSetCameraDirectionAP()
{
  const double vn[3] = { 0, 1, 0 };
  const double vup[3] = { 0, 0, 1 };
  this->SetCamera( vn, vup );
}

void 
QtVtkViewer::slotSetCameraDirectionLR()
{
  const double vn[3] = { 1, 0, 0 };
  const double vup[3] = { 0, 0, 1 };
  this->SetCamera( vn, vup );
}

void 
QtVtkViewer::slotSetCameraDirectionHF()
{
  const double vn[3] = { 0, 0, 1 };
  const double vup[3] = { 0, 1, 0 };
  this->SetCamera( vn, vup );
}

void 
QtVtkViewer::slotSetCameraDirectionPA()
{
  const double vn[3] = { 0, -1, 0 };
  const double vup[3] = { 0, 0, 1 };
  this->SetCamera( vn, vup );
}

void 
QtVtkViewer::slotSetCameraDirectionRL()
{
  const double vn[3] = { -1, 0, 0 };
  const double vup[3] = { 0, 0, 1 };
  this->SetCamera( vn, vup );
}

void 
QtVtkViewer::slotSetCameraDirectionFH()
{
  const double vn[3] = { 0, 0, -1 };
  const double vup[3] = { 0, 1, 0 };
  this->SetCamera( vn, vup );
}

void 
QtVtkViewer::slotAzimuth( int angle )
{
  vtkCamera* camera = Renderer[1][1]->GetActiveCamera();
  camera->Azimuth( angle );
  camera->OrthogonalizeViewUp();
  Renderer[1][1]->ResetCameraClippingRange();
  RenderWindow->Render();
}

void
QtVtkViewer::slotElevation( int angle )
{
  vtkCamera* camera = Renderer[1][1]->GetActiveCamera();
  camera->Elevation( angle );
  Renderer[1][1]->ResetCameraClippingRange();
  RenderWindow->Render();
}

void
QtVtkViewer::slotDolly( float dolly )
{
  for ( int row = 0; row < 2; ++row )
    for ( int col = 0; col < 2; ++col ) 
      {
      vtkCamera* camera = Renderer[row][col]->GetActiveCamera();
      camera->Dolly( dolly );
      Renderer[row][col]->ResetCameraClippingRange();
      }
  RenderWindow->Render();
}

void
QtVtkViewer::SetCamera
( const double* vn, const double* vup, const int row, const int col )
{
  vtkCamera* camera = Renderer[row][col]->GetActiveCamera();

  double bounds[6];
  Renderer[row][col]->ComputeVisiblePropBounds( bounds );

  double center[3];
  double distance;
  double width;

  center[0] = (bounds[0] + bounds[1])/2.0;
  center[1] = (bounds[2] + bounds[3])/2.0;
  center[2] = (bounds[4] + bounds[5])/2.0;

  width = std::max( std::max( bounds[5] - bounds[4], bounds[3] - bounds[2] ), bounds[1] - bounds[0] );
  
  // If we have just a single point, pick a width of 1.0
  width = (width==0)?(1.0):(width);

  // angle = 2*atan((h/2)/d) 
  distance = width / (tan(camera->GetViewAngle() * M_PI/360.0));
  distance = std::max( distance, width/2.0 );

  // check view-up vector against view plane normal
  camera->SetViewUp( vup );

  // update the camera
  camera->SetFocalPoint(center[0],center[1],center[2]);
  camera->SetPosition(center[0]+distance*vn[0], center[1]+distance*vn[1], center[2]+distance*vn[2]);
  camera->OrthogonalizeViewUp();

  Renderer[row][col]->ResetCameraClippingRange();
  RenderWindow->Render();
}
  
void
QtVtkViewer::slotSetLayoutSingle()
{
  for ( int row = 0; row < 2; ++row )
    for ( int col = 0; col < 2; ++col )
      if ( row && col )
	Renderer[row][col]->SetViewport( 0, 0, 1, 1 );
      else
	RenderWindow->RemoveRenderer( Renderer[row][col] );      
  this->slotUpdateRenderers();
}

void
QtVtkViewer::slotSetLayoutPanel()
{
  for ( int row = 0; row < 2; ++row )
    for ( int col = 0; col < 2; ++col ) 
      {
      Renderer[row][col]->SetViewport(  0.5 * col, 0.5 * (1-row), 0.5 * (col+1), 0.5 * (2-row) );
      }
  this->slotUpdateRenderers();
}

void
QtVtkViewer::slotUpdateRenderers()
{
  RenderWindow->Render();
}

void
QtVtkViewer::slotExportImage( const QString& filename )
{
  if ( ! WindowToImage ) 
    {
    WindowToImage = vtkWindowToImageFilter::New();
    WindowToImage->SetInput( RenderWindow );
    }
  
  if ( ! TIFFWriter ) 
    {
    TIFFWriter = vtkTIFFWriter::New();
    TIFFWriter->SetInput( WindowToImage->GetOutput() );
    }

  WindowToImage->Modified();
  TIFFWriter->SetFileName( filename.latin1() );
  TIFFWriter->Write();
}

} // namespace

#endif // #if defined(CMTK_HAVE_QT) && defined(CMTK_HAVE_VTK)
