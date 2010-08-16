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

#include <QtGui/QActionGroup>
#include <QtGui/QPainter>
#include <QtGui/QImage>

cmtk::FusionViewApplication
::FusionViewApplication( int argc, char* argv[] ) 
  : QApplication( argc, argv ),
    m_MainWindow( new QMainWindow ),
    m_SliceAxis( 2 )
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

  this->m_SliceIndex = this->m_FixedVolume->GetDims()[this->m_SliceAxis] / 2;  
  this->m_MainWindowUI.sliceSlider->setRange( 0, this->m_FixedVolume->GetDims()[this->m_SliceAxis] );
  this->m_MainWindowUI.sliceSlider->setValue( this->m_SliceIndex );
  QObject::connect( this->m_MainWindowUI.sliceSlider, SIGNAL( valueChanged( int ) ), this, SLOT( SetFixedSlice( int ) ) );

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
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom25 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom50 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom100 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom200 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom300 );
  zoomGroup->addAction( this->m_MainWindowUI.actionZoom400 );
  
  QActionGroup* sliceGroup = new QActionGroup( this->m_MainWindow );
  sliceGroup->setExclusive( true );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceAxial_XY );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceCoronal_XZ );
  sliceGroup->addAction( this->m_MainWindowUI.actionSliceSagittal_YZ );
  
  QActionGroup* interpGroup = new QActionGroup( this->m_MainWindow );
  interpGroup->setExclusive( true );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpLinear );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpCubic );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpSinc );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpNearestNeighbour );
  interpGroup->addAction( this->m_MainWindowUI.actionInterpPartialVolume );
  
  this->m_MainWindow->show();
}

void
cmtk::FusionViewApplication
::SetFixedSlice( int slice )
{
  if ( this->m_SliceIndex != slice )
    {
    this->m_SliceIndex = slice;
    this->m_FixedSlice = this->m_FixedVolume->ExtractSlice( this->m_SliceAxis, this->m_SliceIndex );
    this->UpdateFixedSlice();
    }
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

  this->UpdateWidget( this->m_MainWindowUI.fixedImgWidget, *(this->m_FixedSlice), this->m_ColorTableFix, this->m_MainWindowUI.blackSliderFix->value(), this->m_MainWindowUI.whiteSliderFix->value() );
}

void
cmtk::FusionViewApplication
::UpdateMovingSlice()
{
}

void
cmtk::FusionViewApplication
::UpdateWidget( QLabel* widget, const UniformVolume& slice, const QVector<QRgb>& colorTable, const float blackLevel, const float whiteLevel )
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
  
  widget->setPixmap( QPixmap::fromImage( image ) );
}
