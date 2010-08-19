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

#include "cmtkFusionViewApplication.h"

#include <System/cmtkCommandLine.h>

#include <IO/cmtkVolumeIO.h>
#include <IO/cmtkXformListIO.h>

#include <Qt/cmtkQtIcons.h>

#include <QtGui/QActionGroup>
#include <QtGui/QPainter>
#include <QtGui/QImage>
#include <QtGui/QScrollBar>

#include <vector>
#include <string>

cmtk::FusionViewApplication
::FusionViewApplication( int argc, char* argv[] ) 
  : QApplication( argc, argv ),
    m_MainWindow( new QMainWindow ),
    m_AffineOnly( false ),
    m_SliceAxis( -1 ),
    m_SliceIndex( -1 ),
    m_Interpolator( Interpolators::LINEAR ),
    m_ZoomFactor( 1.0 ),
    m_Transparency( 1.0 ),
    m_CursorDisplayed( true )
{
  CommandLine cl;
  cl.SetProgramInfo( CommandLine::PRG_TITLE, "Fusion viewer." );

  const char* imagePathFix;
  cl.AddParameter( &imagePathFix, "FixedImage", "Fixed image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

  const char* imagePathMov;
  cl.AddParameter( &imagePathMov, "MovingImage", "Moving image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

  std::vector<std::string> xformList;
  cl.AddParameterVector( &xformList, "XformList", "List of concatenated transformations. Insert '--inverse' to use the inverse of the transformation listed next." )
    ->SetProperties( cmtk::CommandLine::PROPS_XFORM | cmtk::CommandLine::PROPS_OPTIONAL );  
  
  try
    {
    cl.Parse( argc, const_cast<const char**>( argv ) );
    }
  catch ( const CommandLine::Exception& ex )
    {
    throw(ex);
    }

  this->m_XformList = XformListIO::MakeFromStringList( xformList );
  this->m_XformListAllAffine = this->m_XformList.MakeAllAffine();

  this->m_Fixed.m_Volume = VolumeIO::ReadOriented( imagePathFix );
  if ( ! this->m_Fixed.m_Volume )
    {
    exit( 1 );
    }
  this->m_Fixed.m_DataRange = this->m_Fixed.m_Volume->GetData()->GetRange();
  this->m_CursorPosition[0] = this->m_Fixed.m_Volume->GetDims()[0] / 2;
  this->m_CursorPosition[1] = this->m_Fixed.m_Volume->GetDims()[1] / 2;

  this->m_Moving.m_Volume = VolumeIO::ReadOriented( imagePathMov );
  if ( ! this->m_Moving.m_Volume )
    {
    exit( 1 );
    }
  this->m_Moving.m_DataRange = this->m_Moving.m_Volume->GetData()->GetRange();

  this->m_MainWindowUI.setupUi( this->m_MainWindow );
  this->m_MainWindow->setWindowIcon( QtIcons::WindowIcon() );

  this->InitViewData( this->m_Fixed, this->m_MainWindowUI.fixedView );
  this->InitViewData( this->m_Moving, this->m_MainWindowUI.movingView );

  QObject::connect( this->m_MainWindowUI.blackSliderFix, SIGNAL( valueChanged( int ) ), this, SLOT( fixedBlackWhiteChanged() ) );
  QObject::connect( this->m_MainWindowUI.whiteSliderFix, SIGNAL( valueChanged( int ) ), this, SLOT( fixedBlackWhiteChanged() ) );
  
  QObject::connect( this->m_MainWindowUI.blackSliderMov, SIGNAL( valueChanged( int ) ), this, SLOT( movingBlackWhiteChanged() ) );
  QObject::connect( this->m_MainWindowUI.whiteSliderMov, SIGNAL( valueChanged( int ) ), this, SLOT( movingBlackWhiteChanged() ) );
  
  QActionGroup* zoomGroup = new QActionGroup( this->m_MainWindow );
  zoomGroup->setExclusive( true );
  this->m_MainWindowUI.actionZoom25->setData( QVariant( 0.25 ) );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom25 );
  this->m_MainWindowUI.actionZoom50->setData( QVariant( 0.50 ) );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom50 );
  this->m_MainWindowUI.actionZoom100->setData( QVariant( 1.00 ) );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom100 );
  this->m_MainWindowUI.actionZoom200->setData( QVariant( 2.00 ) );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom200 );
  this->m_MainWindowUI.actionZoom300->setData( QVariant( 3.00 ) );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom300 );
  this->m_MainWindowUI.actionZoom400->setData( QVariant( 4.00 ) );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom400 );

  QObject::connect( zoomGroup, SIGNAL( triggered( QAction* ) ), this, SLOT( changeZoom( QAction* ) ) );
  
  QActionGroup* sliceGroup = new QActionGroup( this->m_MainWindow );
  sliceGroup->setExclusive( true );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceAxial_XY );
  this->m_MainWindowUI.actionSliceAxial_XY->setData( QVariant( AXIS_Z ) );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceCoronal_XZ );
  this->m_MainWindowUI.actionSliceCoronal_XZ->setData( QVariant( AXIS_Y ) );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceSagittal_YZ );
  this->m_MainWindowUI.actionSliceSagittal_YZ->setData( QVariant( AXIS_X ) );
  
  QObject::connect( sliceGroup, SIGNAL( triggered( QAction* ) ), this, SLOT( changeSliceDirection( QAction* ) ) );
  
  QActionGroup* interpGroup = new QActionGroup( this->m_MainWindow );
  interpGroup->setExclusive( true );
  this->m_MainWindowUI.actionInterpLinear->setData( QVariant( Interpolators::LINEAR ) );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpLinear );
  this->m_MainWindowUI.actionInterpCubic->setData( QVariant( Interpolators::CUBIC ) );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpCubic );
  this->m_MainWindowUI.actionInterpSinc->setData( QVariant( Interpolators::COSINE_SINC ) );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpSinc );
  this->m_MainWindowUI.actionInterpNearestNeighbour->setData( QVariant( Interpolators::NEAREST_NEIGHBOR ) );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpNearestNeighbour );
  this->m_MainWindowUI.actionInterpPartialVolume->setData( QVariant( Interpolators::PARTIALVOLUME ) );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpPartialVolume );

  QObject::connect( interpGroup, SIGNAL( triggered( QAction* ) ), this, SLOT( changeInterpolator( QAction* ) ) );
  
  this->m_MainWindowUI.alphaSlider->setRange( 0, 1000 );
  this->m_MainWindowUI.alphaSlider->setValue( this->m_Transparency * 1000 );
  QObject::connect( this->m_MainWindowUI.alphaSlider, SIGNAL( valueChanged( int ) ), this, SLOT( setTransparency( int ) ) );
  
  this->changeSliceDirection( AXIS_Z );
  QObject::connect( this->m_MainWindowUI.sliceSlider, SIGNAL( valueChanged( int ) ), this, SLOT( setFixedSlice( int ) ) );

  QObject::connect( this->m_MainWindowUI.actionLinkedCursor, SIGNAL( toggled( bool ) ), this, SLOT( setLinkedCursorFlag( bool ) ) );  

  QObject::connect( this->m_MainWindowUI.actionAffineOnly, SIGNAL( toggled( bool ) ), this, SLOT( setAffineOnly( bool ) ) );  
    
  // synchronize sliders of the two graphics views
  QObject::connect( this->m_MainWindowUI.fixedView->horizontalScrollBar(), SIGNAL( valueChanged( int ) ), 
		    this->m_MainWindowUI.movingView->horizontalScrollBar(), SLOT( setValue( int ) ) );
  QObject::connect( this->m_MainWindowUI.fixedView->verticalScrollBar(), SIGNAL( valueChanged( int ) ), 
		    this->m_MainWindowUI.movingView->verticalScrollBar(), SLOT( setValue( int) ) );

  QObject::connect( this->m_MainWindowUI.movingView->horizontalScrollBar(), SIGNAL( valueChanged( int ) ), 
		    this->m_MainWindowUI.fixedView->horizontalScrollBar(), SLOT( setValue( int ) ) );
  QObject::connect( this->m_MainWindowUI.movingView->verticalScrollBar(), SIGNAL( valueChanged( int ) ), 
		    this->m_MainWindowUI.fixedView->verticalScrollBar(), SLOT( setValue( int ) ) );

  this->m_MainWindow->show();
}

void
cmtk::FusionViewApplication
::InitViewData( Self::Data& data, QGraphicsView* view )
{
  const QPen redPen( QColor( 255, 0, 0 ) );

  data.m_Scene = new QGraphicsScene;
  (data.m_CursorLines[0] = data.m_Scene->addLine( 0, 0, 0, 0, redPen ))->setZValue( 100 );
  (data.m_CursorLines[1] = data.m_Scene->addLine( 0, 0, 0, 0, redPen ))->setZValue( 100 );
  (data.m_PixmapItem = new QGraphicsPixmapItemEvents)->setZValue( 0 );
  data.m_Scene->addItem( data.m_PixmapItem );

  QObject::connect( data.m_PixmapItem, SIGNAL( mousePressed( QGraphicsSceneMouseEvent* ) ), this, SLOT( mousePressed( QGraphicsSceneMouseEvent* ) ) );

  data.m_View = view;
  data.m_View->setScene( data.m_Scene );
}


void
cmtk::FusionViewApplication
::setTransparency( int slice )
{
  this->m_Transparency = static_cast<float>( slice ) / 1000;
  this->UpdateMovingImage();
}

void
cmtk::FusionViewApplication
::setLinkedCursorFlag( bool flag )
{
  this->m_CursorDisplayed = flag;

  for ( int i = 0; i < 2; ++i )
    {
    this->m_Fixed.m_CursorLines[i]->setVisible( this->m_CursorDisplayed );
    this->m_Moving.m_CursorLines[i]->setVisible( this->m_CursorDisplayed );
    }
}

void
cmtk::FusionViewApplication
::setFixedSlice( int slice )
{
  if ( this->m_SliceIndex != slice )
    {
    this->m_SliceIndex = slice;
    this->m_Fixed.m_Slice = this->m_Fixed.m_Volume->ExtractSlice( this->m_SliceAxis, this->m_SliceIndex );
    this->UpdateFixedImage();

    this->UpdateMovingSlice();

    this->m_MainWindowUI.sliceLabel->setText( QString("Slice: %1").arg( this->m_SliceIndex ) );
    }
}

void
cmtk::FusionViewApplication
::changeZoom( QAction* action )
{
  this->m_ZoomFactor = static_cast<float>( action->data().toDouble() ); // older Qt doesn't have QVariant::toFloat()
  this->UpdateFixedImage();
  this->UpdateMovingImage();
}

void
cmtk::FusionViewApplication
::changeInterpolator( QAction* action )
{
  this->m_Interpolator = static_cast<Interpolators::InterpolationEnum>( action->data().toInt() );
  this->UpdateMovingSlice();
}

void
cmtk::FusionViewApplication
::setAffineOnly( bool affineOnly )
{
  this->m_AffineOnly = affineOnly;
  this->UpdateMovingSlice();
}

void
cmtk::FusionViewApplication
::fixedBlackWhiteChanged()
{
  this->UpdateFixedImage();
}

void
cmtk::FusionViewApplication
::movingBlackWhiteChanged()
{
  this->UpdateMovingImage();
}

void
cmtk::FusionViewApplication
::changeSliceDirection( QAction* action )
{
  this->changeSliceDirection( action->data().toInt() );
}

void
cmtk::FusionViewApplication
::mousePressed( QGraphicsSceneMouseEvent* event )  
{
  this->m_CursorPosition[0] = event->pos().x();
  this->m_CursorPosition[1] = event->pos().y();

  if ( this->m_CursorDisplayed )
    {
    this->UpdateView( this->m_Fixed, this->m_Fixed.m_Image );
    this->UpdateView( this->m_Moving, this->m_FusedImage );
    }
}

void
cmtk::FusionViewApplication
::changeSliceDirection( const int sliceAxis )
{
  if ( sliceAxis != this->m_SliceAxis )
    {
    this->m_SliceAxis = sliceAxis;

    this->m_SliceIndex = -1; // unset previously set slice index to ensure update of moving slice
    this->m_MainWindowUI.sliceSlider->setRange( 0, this->m_Fixed.m_Volume->GetDims()[this->m_SliceAxis]-1 );
    this->setFixedSlice( this->m_Fixed.m_Volume->GetDims()[this->m_SliceAxis] / 2 ); 
    this->m_MainWindowUI.sliceSlider->setValue( this->m_SliceIndex );

    this->UpdateMovingSlice();

    const char* labelFrom[3] = { "Left", "Posterior", "Inferior" };
    const char* labelTo[3] = { "Right", "Anterior", "Superior" };

    this->m_MainWindowUI.sliceLabelFrom->setText( labelFrom[this->m_SliceAxis] );
    this->m_MainWindowUI.sliceLabelTo->setText( labelTo[this->m_SliceAxis] );
    }
}

void
cmtk::FusionViewApplication
::UpdateMovingSlice()
{
  ReformatVolume::Plain plain( TYPE_FLOAT );
  UniformVolumeInterpolatorBase::SmartPtr interpolator ( ReformatVolume::CreateInterpolator( this->m_Interpolator, this->m_Moving.m_Volume ) );
  
  const XformList noXforms;
  TypedArray::SmartPtr reformatData( ReformatVolume::ReformatUnmasked( this->m_Fixed.m_Slice, this->m_AffineOnly ? this->m_XformListAllAffine :  this->m_XformList, noXforms, plain, this->m_Moving.m_Volume, interpolator ) );
  
  UniformVolume::SmartPtr movingSlice = this->m_Fixed.m_Slice->CloneGrid();
  movingSlice->SetData( reformatData );
  
  this->m_Moving.m_Slice = movingSlice;
  this->UpdateMovingImage();
}

void
cmtk::FusionViewApplication
::UpdateFixedImage()
{
  this->m_Fixed.m_ColorTable.resize( 256 );
  for ( int i = 0; i < 256; ++i )
    {
    this->m_Fixed.m_ColorTable[i] = QColor( i, i, i ).rgb();
    }

  const float black = this->m_Fixed.m_DataRange.m_LowerBound + this->m_Fixed.m_DataRange.Width() * this->m_MainWindowUI.blackSliderFix->value() / 500;
  const float white = this->m_Fixed.m_DataRange.m_LowerBound + this->m_Fixed.m_DataRange.Width() * this->m_MainWindowUI.whiteSliderFix->value() / 500;

  this->MakeImage( this->m_Fixed.m_Image, *(this->m_Fixed.m_Slice), this->m_Fixed.m_ColorTable, black, white );
  this->UpdateView( this->m_Fixed, this->m_Fixed.m_Image );
}

void
cmtk::FusionViewApplication
::UpdateMovingImage()
{
  this->m_Moving.m_ColorTable.resize( 256 );
  for ( int i = 0; i < 256; ++i )
    {
    this->m_Moving.m_ColorTable[i] = QColor( i, i, i ).rgb();
    }

  const float black = this->m_Moving.m_DataRange.m_LowerBound + this->m_Moving.m_DataRange.Width() * this->m_MainWindowUI.blackSliderMov->value() / 500;
  const float white = this->m_Moving.m_DataRange.m_LowerBound + this->m_Moving.m_DataRange.Width() * this->m_MainWindowUI.whiteSliderMov->value() / 500;

  this->MakeImage( this->m_Moving.m_Image, *(this->m_Moving.m_Slice), this->m_Moving.m_ColorTable, black, white );

  this->m_FusedImage = QImage( this->m_Moving.m_Image.width(), this->m_Moving.m_Image.height(), QImage::Format_RGB32 );
  for ( int y = 0; y < this->m_Moving.m_Image.height(); ++y )
    {
    for ( int x = 0; x < this->m_Moving.m_Image.width(); ++x )
      {
      QColor rgbMov( this->m_Moving.m_Image.pixel( x, y ) );
      const QColor rgbFix( this->m_Fixed.m_Image.pixel( x, y ) );

      rgbMov = QColor( this->m_Transparency * rgbMov.red() + (1.0-this->m_Transparency) * rgbFix.red(),
		       this->m_Transparency * rgbMov.green() + (1.0-this->m_Transparency) * rgbFix.green(),
		       this->m_Transparency * rgbMov.blue() + (1.0-this->m_Transparency) * rgbFix.blue() );
      this->m_FusedImage.setPixel( x, y, rgbMov.rgb() );
      }
    }

  this->UpdateView( this->m_Moving, this->m_FusedImage );
}

void
cmtk::FusionViewApplication
::MakeImage( QImage& image, const UniformVolume& slice, const QVector<QRgb>& colorTable, const float blackLevel, const float whiteLevel )
{
  // table: which are the x and y directions for slices along the three orthogonal orientations?
  const int idxXtable[3] = { 1, 0, 0 };
  const int idxYtable[3] = { 2, 2, 1 };

  const int idxX = idxXtable[this->m_SliceAxis];
  const int idxY = idxYtable[this->m_SliceAxis];
  
  int dimX = slice.GetDims()[idxX];
  int dimY = slice.GetDims()[idxY];

  const float dX = slice.Deltas()[idxX];
  const float dY = slice.Deltas()[idxY];

  // make sure pixels are displayed square
  this->m_ScalePixels[0] = std::max<float>( 1.0, dX / dY );
  this->m_ScalePixels[1] = std::max<float>( 1.0, dY / dX );
  
  image = QImage( dimX, dimY, QImage::Format_Indexed8 );
  image.setColorTable( colorTable );
  
  const float scaleLevel = 1.0 / (whiteLevel-blackLevel);
  
  size_t idx = 0;
  for ( int y = 0; y < dimY; ++y )
    {
    for ( int x = 0; x < dimX; ++x, ++idx )
      {
      image.setPixel( x, y, static_cast<int>( 255 * std::min<float>( 1, std::max<float>( 0, (slice.GetDataAt( idx ) - blackLevel) * scaleLevel ) ) ) );
      }
    }
}

void
cmtk::FusionViewApplication
::UpdateView( Self::Data& data, QImage& image )
{
  data.m_PixmapItem->setPixmap( QPixmap::fromImage( image ) );

  const QRectF bb = data.m_PixmapItem->boundingRect();
  data.m_Scene->setSceneRect( bb );

  if ( this->m_CursorDisplayed )
    {
    data.m_CursorLines[0]->setLine( this->m_CursorPosition[0], bb.top(), this->m_CursorPosition[0], bb.bottom() );
    data.m_CursorLines[1]->setLine( bb.left(), this->m_CursorPosition[1], bb.right(), this->m_CursorPosition[1] );
    }

  QTransform zoomTransform = QTransform::fromScale( -this->m_ZoomFactor * this->m_ScalePixels[0], -this->m_ZoomFactor * this->m_ScalePixels[1] );
  data.m_View->setTransform( zoomTransform );
}
