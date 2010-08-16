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
#include <IO/cmtkXformIO.h>

#include <QtGui/QActionGroup>
#include <QtGui/QPainter>
#include <QtGui/QImage>

cmtk::FusionViewApplication
::FusionViewApplication( int argc, char* argv[] ) 
  : QApplication( argc, argv ),
    m_MainWindow( new QMainWindow ),
    m_SliceAxis( 2 ),
    m_Interpolator( Interpolators::LINEAR ),
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

    const char* next = cl.GetNextOptional();
    while ( next ) 
      {
      if ( ! strcmp( next, "-j" ) || ! strcmp( next, "--jacobian" ) )
	break;
      
      const bool inverse = ! strcmp( next, "-i" ) || ! strcmp( next, "--inverse" );
      if ( inverse ) 
	next = cl.GetNext();
      
      Xform::SmartPtr xform( XformIO::Read( next ) );
      if ( ! xform ) 
	{
	cmtk::StdErr << "ERROR: could not read target-to-reference transformation from " << next << "\n";
	exit( 1 );
	}
      
      this->m_XformList.Add( xform, inverse );
      next = cl.GetNextOptional();
      }
    }
  catch ( const CommandLine::Exception& ex )
    {
    throw(ex);
    }

  this->m_FixedVolume = VolumeIO::Read( imagePathFix );
  if ( ! this->m_FixedVolume )
    {
    exit( 1 );
    }

  this->m_MovingVolume = VolumeIO::Read( imagePathMov );
  if ( ! this->m_MovingVolume )
    {
    exit( 1 );
    }

  this->m_MainWindowUI.setupUi( this->m_MainWindow );

  this->m_MainWindowUI.alphaSlider->setRange( 0, 1000 );

  const Types::DataItemRange rangeFix = this->m_FixedVolume->GetData()->GetRange();
  this->m_MainWindowUI.blackSliderFix->setRange( rangeFix.m_LowerBound, rangeFix.m_UpperBound );
  this->m_MainWindowUI.blackSliderFix->setValue( rangeFix.m_LowerBound );
  this->m_MainWindowUI.whiteSliderFix->setRange( rangeFix.m_LowerBound, rangeFix.m_UpperBound );
  this->m_MainWindowUI.whiteSliderFix->setValue( rangeFix.m_UpperBound );
  
  const Types::DataItemRange rangeMov = this->m_MovingVolume->GetData()->GetRange();
  this->m_MainWindowUI.blackSliderMov->setRange( rangeMov.m_LowerBound, rangeMov.m_UpperBound );
  this->m_MainWindowUI.blackSliderMov->setValue( rangeMov.m_LowerBound );
  this->m_MainWindowUI.whiteSliderMov->setRange( rangeMov.m_LowerBound, rangeMov.m_UpperBound );
  this->m_MainWindowUI.whiteSliderMov->setValue( rangeMov.m_UpperBound );
  
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
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceCoronal_XZ );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceSagittal_YZ );
  
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
  
  this->m_SliceIndex = this->m_FixedVolume->GetDims()[this->m_SliceAxis] / 2;  
  this->m_MainWindowUI.sliceSlider->setRange( 0, this->m_FixedVolume->GetDims()[this->m_SliceAxis]-1 );
  this->m_MainWindowUI.sliceSlider->setValue( this->m_SliceIndex );
  QObject::connect( this->m_MainWindowUI.sliceSlider, SIGNAL( valueChanged( int ) ), this, SLOT( setFixedSlice( int ) ) );

  this->setFixedSlice( this->m_MainWindowUI.sliceSlider->value() );

  this->m_MainWindow->show();
}

void
cmtk::FusionViewApplication
::setFixedSlice( int slice )
{
  if ( this->m_SliceIndex != slice )
    {
    this->m_SliceIndex = slice;
    this->m_FixedSlice = this->m_FixedVolume->ExtractSlice( this->m_SliceAxis, this->m_SliceIndex );
    this->UpdateFixedSlice();

    ReformatVolume::Plain plain( TYPE_FLOAT );
    UniformVolumeInterpolatorBase::SmartPtr interpolator ( ReformatVolume::CreateInterpolator( this->m_Interpolator, this->m_MovingVolume ) );

    const XformList noXforms;
    TypedArray::SmartPtr reformatData( ReformatVolume::Reformat( this->m_FixedSlice, this->m_XformList, noXforms, plain, this->m_MovingVolume, interpolator ) );

    UniformVolume::SmartPtr movingSlice = this->m_FixedSlice->CloneGrid();
    movingSlice->SetData( reformatData );

    this->m_MovingSlice = movingSlice;
    this->UpdateMovingSlice();
    }
}

void
cmtk::FusionViewApplication
::changeZoom( QAction* action )
{
  this->m_ZoomFactor = action->data().toFloat();
  this->UpdateFixedSlice();
  this->UpdateMovingSlice();
}

void
cmtk::FusionViewApplication
::changeInterpolator( QAction* action )
{
  this->m_Interpolator = static_cast<Interpolators::InterpolationEnum>( action->data().toInt() );
  this->UpdateFixedSlice();
  this->UpdateMovingSlice();
}

void
cmtk::FusionViewApplication
::UpdateFixedSlice()
{
  this->m_ColorTableFix.resize( 256 );
  for ( int i = 0; i < 256; ++i )
    {
    this->m_ColorTableFix[i] = QColor( i, i, i ).rgb();
    }

  this->UpdateView( this->m_MainWindowUI.fixedView, *(this->m_FixedSlice), this->m_ColorTableFix, this->m_MainWindowUI.blackSliderFix->value(), this->m_MainWindowUI.whiteSliderFix->value() );
}

void
cmtk::FusionViewApplication
::UpdateMovingSlice()
{
  this->m_ColorTableMov.resize( 256 );
  for ( int i = 0; i < 256; ++i )
    {
    this->m_ColorTableMov[i] = QColor( i, i, i ).rgb();
    }

  this->UpdateView( this->m_MainWindowUI.movingView, *(this->m_MovingSlice), this->m_ColorTableMov, this->m_MainWindowUI.blackSliderMov->value(), this->m_MainWindowUI.whiteSliderMov->value() );
}

void
cmtk::FusionViewApplication
::UpdateView( QGraphicsView* view, const UniformVolume& slice, const QVector<QRgb>& colorTable, const float blackLevel, const float whiteLevel )
{
  QImage image( slice.GetDims()[0], slice.GetDims()[1], QImage::Format_Indexed8 );
  image.setColorTable( colorTable );

  const float scaleLevel = 1.0 / (whiteLevel-blackLevel);
  
  size_t idx = 0;
  for ( int y = 0; y < slice.GetDims()[1]; ++y )
    {
    for ( int x = 0; x < slice.GetDims()[0]; ++x, ++idx )
      {
      image.setPixel( x, y, static_cast<int>( 255 * std::min<float>( 1, std::max<float>( 0, (slice.GetDataAt( idx ) - blackLevel) * scaleLevel ) ) ) );
      }
    }
  
  QGraphicsScene* scene = new QGraphicsScene;
  scene->addPixmap( QPixmap::fromImage( image ) );

  QTransform zoomTransform = QTransform::fromScale( this->m_ZoomFactor, -this->m_ZoomFactor );
  view->setTransform( zoomTransform );

  view->setScene( scene );
}
