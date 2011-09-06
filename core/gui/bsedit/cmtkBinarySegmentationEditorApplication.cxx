/*
//
//  Copyright 2010-2011 SRI International
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

#include "cmtkBinarySegmentationEditorApplication.h"

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

cmtk::BinarySegmentationEditorApplication
::BinarySegmentationEditorApplication( int& argc, char* argv[] ) 
  : QApplication( argc, argv ),
    m_MainWindow( new QMainWindow ),
    m_SliceAxis( -1 ),
    m_SliceIndex( -1 ),
    m_ZoomFactor( 1.0 )
{
  CommandLine cl;
  cl.SetProgramInfo( CommandLine::PRG_TITLE, "Fusion viewer." );

  const char* imagePathFix;
  cl.AddParameter( &imagePathFix, "FixedImage", "Fixed image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );
  
  const char* imagePathMov;
  cl.AddParameter( &imagePathMov, "MovingImage", "Moving image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

  try
    {
    cl.Parse( argc, const_cast<const char**>( argv ) );
    }
  catch ( const CommandLine::Exception& ex )
    {
    throw(ex);
    }

  this->m_Fixed.m_Volume = VolumeIO::ReadOriented( imagePathFix );
  if ( ! this->m_Fixed.m_Volume )
    {
    StdErr << "Fixed image '" << imagePathFix << "' could not be read.\n";
    throw( ExitException( 1 ) );
    }
  this->m_Fixed.m_DataRange = this->m_Fixed.m_Volume->GetData()->GetRange();

  // per-dimension scale factors make sure non-square pixels are displayed square
  this->m_ScalePixels = this->m_Fixed.m_Volume->Deltas();
  this->m_ScalePixels *= 1.0 / this->m_ScalePixels.MinValue();

  this->m_MainWindowUI.setupUi( this->m_MainWindow );
  this->m_MainWindow->setWindowIcon( QtIcons::WindowIcon() );

  this->InitViewData( this->m_Fixed, this->m_MainWindowUI.fixedView );

  QObject::connect( this->m_MainWindowUI.blackSliderFix, SIGNAL( valueChanged( int ) ), this, SLOT( fixedBlackWhiteChanged() ) );
  QObject::connect( this->m_MainWindowUI.whiteSliderFix, SIGNAL( valueChanged( int ) ), this, SLOT( fixedBlackWhiteChanged() ) );
  
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

  QActionGroup* fixedColorGroup = new QActionGroup( this->m_MainWindow );
  fixedColorGroup->setExclusive( true );
  this->m_MainWindowUI.actionFixedGrey->setData( QVariant( 0 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedGrey );
  this->m_MainWindowUI.actionFixedRed->setData( QVariant( 1 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedRed );
  this->m_MainWindowUI.actionFixedGreen->setData( QVariant( 2 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedGreen );
  this->m_MainWindowUI.actionFixedBlue->setData( QVariant( 3 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedBlue );
  this->m_MainWindowUI.actionFixedCyan->setData( QVariant( 4 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedCyan );
  this->m_MainWindowUI.actionFixedYellow->setData( QVariant( 5 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedYellow );
  this->m_MainWindowUI.actionFixedMagenta->setData( QVariant( 6 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedMagenta );
  this->m_MainWindowUI.actionFixedBlueRed->setData( QVariant( 7 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedBlueRed );
  this->m_MainWindowUI.actionFixedRedBlue->setData( QVariant( 8 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedRedBlue );
  this->m_MainWindowUI.actionFixedLabels->setData( QVariant( 9 ) );
  fixedColorGroup->addAction( this->m_MainWindowUI.actionFixedLabels );

  QObject::connect( fixedColorGroup, SIGNAL( triggered( QAction* ) ), this, SLOT( changeFixedColor( QAction* ) ) );

  QActionGroup* sliceGroup = new QActionGroup( this->m_MainWindow );
  sliceGroup->setExclusive( true );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceAxial_XY );
  this->m_MainWindowUI.actionSliceAxial_XY->setData( QVariant( AXIS_Z ) );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceCoronal_XZ );
  this->m_MainWindowUI.actionSliceCoronal_XZ->setData( QVariant( AXIS_Y ) );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceSagittal_YZ );
  this->m_MainWindowUI.actionSliceSagittal_YZ->setData( QVariant( AXIS_X ) );
  
  QObject::connect( sliceGroup, SIGNAL( triggered( QAction* ) ), this, SLOT( changeSliceDirection( QAction* ) ) );
  
  this->changeSliceDirection( AXIS_Z );
  QObject::connect( this->m_MainWindowUI.sliceSlider, SIGNAL( valueChanged( int ) ), this, SLOT( setFixedSlice( int ) ) );

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
cmtk::BinarySegmentationEditorApplication
::InitViewData( Self::Data& data, QGraphicsView* view )
{
  data.m_ColorMapIndex = 0;

  const QPen redPen( QColor( 255, 0, 0 ) );

  data.m_Scene = new QGraphicsScene;
  QObject::connect( data.m_PixmapItem, SIGNAL( mousePressed( QGraphicsSceneMouseEvent* ) ), this, SLOT( mousePressed( QGraphicsSceneMouseEvent* ) ) );

  data.m_View = view;
  data.m_View->setScene( data.m_Scene );
}


void
cmtk::BinarySegmentationEditorApplication
::changeFixedColor( QAction* action )
{
  this->m_Fixed.m_ColorMapIndex = action->data().toInt();
  this->UpdateFixedImage();
}

void
cmtk::BinarySegmentationEditorApplication
::setFixedSlice( int slice )
{
  if ( this->m_SliceIndex != slice )
    {
    this->m_SliceIndex = slice;

    this->m_Fixed.m_Slice = this->m_Fixed.m_Volume->ExtractSlice( this->m_SliceAxis, this->m_SliceIndex );
    this->UpdateFixedImage();

    this->m_MainWindowUI.sliceLabel->setText( QString("Slice: %1").arg( this->m_SliceIndex ) );
    }
}

void
cmtk::BinarySegmentationEditorApplication
::changeZoom( QAction* action )
{
  this->m_ZoomFactor = static_cast<float>( action->data().toDouble() ); // older Qt doesn't have QVariant::toFloat()
  this->UpdateFixedImage();
}

void
cmtk::BinarySegmentationEditorApplication
::fixedBlackWhiteChanged()
{
  this->UpdateFixedImage();
}

void
cmtk::BinarySegmentationEditorApplication
::changeSliceDirection( QAction* action )
{
  this->changeSliceDirection( action->data().toInt() );
}

void
cmtk::BinarySegmentationEditorApplication
::mousePressed( QGraphicsSceneMouseEvent* event )  
{
  const int idxX = this->GetAxis2DX();
  const int idxY = this->GetAxis2DY();

}

void
cmtk::BinarySegmentationEditorApplication
::changeSliceDirection( const int sliceAxis )
{
  if ( sliceAxis != this->m_SliceAxis )
    {
    this->m_SliceAxis = sliceAxis;

    this->m_SliceIndex = -1; // unset previously set slice index to ensure update of moving slice
    const int newSliceIndex = static_cast<int>( this->m_Fixed.m_Volume->GetDims()[this->m_SliceAxis] / 2 );
    this->m_MainWindowUI.sliceSlider->setRange( 0, this->m_Fixed.m_Volume->GetDims()[this->m_SliceAxis]-1 );
    
    this->setFixedSlice( newSliceIndex );
    this->m_MainWindowUI.sliceSlider->setValue( this->m_SliceIndex );
    
    const char* labelFrom[3] = { "Left", "Posterior", "Inferior" };
    const char* labelTo[3] = { "Right", "Anterior", "Superior" };

    this->m_MainWindowUI.sliceLabelFrom->setText( labelFrom[this->m_SliceAxis] );
    this->m_MainWindowUI.sliceLabelTo->setText( labelTo[this->m_SliceAxis] );
    }
}

void
cmtk::BinarySegmentationEditorApplication
::MakeColorTable( Self::Data& data )
{
  data.m_ColorTable.resize( 256 );

  switch ( data.m_ColorMapIndex )
    {
    default:
    case 0: // black/white
      for ( int i = 0; i < 256; ++i )
	{
	data.m_ColorTable[i] = QColor( i, i, i ).rgb();
	}
      break;
    case 1: // red
      for ( int i = 0; i < 256; ++i )
	{
	data.m_ColorTable[i] = QColor( i, 0, 0 ).rgb();
	}
      break;
    case 2: // green
      for ( int i = 0; i < 256; ++i )
	{
	data.m_ColorTable[i] = QColor( 0, i, 0 ).rgb();
	}
      break;
    case 3: // blue
      for ( int i = 0; i < 256; ++i )
	{
	data.m_ColorTable[i] = QColor( 0, 0, i ).rgb();
	}
      break;
    case 4: // cyan
      for ( int i = 0; i < 256; ++i )
	{
	data.m_ColorTable[i] = QColor( 0, i, i ).rgb();
	}
      break;
    case 5: // yellow
      for ( int i = 0; i < 256; ++i )
	{
	data.m_ColorTable[i] = QColor( i, i, 0 ).rgb();
	}
      break;
    case 6: // magenta
      for ( int i = 0; i < 256; ++i )
	{
	data.m_ColorTable[i] = QColor( i, 0, i ).rgb();
	}
      break;
    case 7: // blue to red
      for ( int i = 0; i < 256; ++i )
	{
	QColor color;
	color.setHsv( 255-i, 255, 255 );
	data.m_ColorTable[i] = color.rgb();
	}
      break;
    case 8: // red to blue
      for ( int i = 0; i < 256; ++i )
	{
	QColor color;
	color.setHsv( i, 255, 255 );
	data.m_ColorTable[i] = color.rgb();
	}
      break;
    }
}

void
cmtk::BinarySegmentationEditorApplication
::UpdateFixedImage()
{
  this->MakeColorTable( this->m_Fixed );

  const float black = this->m_Fixed.m_DataRange.m_LowerBound + this->m_Fixed.m_DataRange.Width() * this->m_MainWindowUI.blackSliderFix->value() / 500;
  const float white = this->m_Fixed.m_DataRange.m_LowerBound + this->m_Fixed.m_DataRange.Width() * this->m_MainWindowUI.whiteSliderFix->value() / 500;

  this->MakeImage( this->m_Fixed.m_Image, *(this->m_Fixed.m_Slice), this->m_Fixed.m_ColorTable, black, white );
  this->UpdateView( this->m_Fixed, this->m_Fixed.m_Image );
}

void
cmtk::BinarySegmentationEditorApplication
::MakeImage( QImage& image, const UniformVolume& slice, const QVector<QRgb>& colorTable, const float blackLevel, const float whiteLevel )
{
  const int idxX = this->GetAxis2DX();
  const int idxY = this->GetAxis2DY();
  
  int dimX = slice.GetDims()[idxX];
  int dimY = slice.GetDims()[idxY];

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
cmtk::BinarySegmentationEditorApplication
::UpdateView( Self::Data& data, QImage& image )
{
  data.m_PixmapItem->setPixmap( QPixmap::fromImage( image ) );

  const QRectF bb = data.m_PixmapItem->boundingRect();
  data.m_Scene->setSceneRect( bb );

  const int idxX = this->GetAxis2DX();
  const int idxY = this->GetAxis2DY();
  
  QTransform zoomTransform = QTransform::fromScale( -this->m_ZoomFactor * this->m_ScalePixels[idxX], -this->m_ZoomFactor * this->m_ScalePixels[idxY] );
  data.m_View->setTransform( zoomTransform );
}
