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

#include <cmtkQtFusionAlpha.h>

#include <qlayout.h>
#include <qmessagebox.h>
#include <QPixmap>

namespace
cmtk
{

QtFusionAlpha::QtFusionAlpha
( QtSimpleFusionApp *const fusionApp, QWidget *const parent, Qt::WFlags flags )
  : 
  QtFusionWindowTemplate( fusionApp, parent, "FusionAlphaWindow", flags ),
  TopIsAlpha( true ),
  StudyTop( NULL ), StudyBottom( NULL )
{
  this->setWindowTitle( "Alpha Blending" );

  StudyNamesBoxTop = new QtStudyNamesBox( this, "StudyBoxTop" );
  StudyNamesBoxTop->slotSetLabel( "Alpha Study:" );
  QObject::connect( StudyNamesBoxTop, SIGNAL( signalActivated(const QString&)), this, SLOT( slotSwitchStudyTop( const QString& ) ) );
  ViewLayout->addWidget( StudyNamesBoxTop );

  StudyNamesBoxBottom = new QtStudyNamesBox( this, "StudyBoxBottom" );
  StudyNamesBoxBottom->slotSetLabel( "Background Study:" );
  QObject::connect( StudyNamesBoxBottom, SIGNAL( signalActivated(const QString&)), this, SLOT( slotSwitchStudyBottom( const QString& ) ) );
  ViewLayout->addWidget( StudyNamesBoxBottom );

  StudyNamesBoxAlpha = new QtStudyNamesBox( this, "StudyBoxAlpha" );
  StudyNamesBoxAlpha->slotSetLabel( "Transparency Study:" );
  QObject::connect( StudyNamesBoxAlpha, SIGNAL( signalActivated(const QString&)), this, SLOT( slotSwitchStudyAlpha( const QString& ) ) );
  ViewLayout->addWidget( StudyNamesBoxAlpha );

  View = new QtScrollRenderView( this );
  ViewLayout->addWidget( View );

  ColormapTop = Colormap::New();
  ColormapTop->SetStandardColormap( PALETTE_GRAY );
  ColormapTop->SetDataRange( 0, 2000 );
  ImageToImageRGBTop = ImageToImageRGB::New();
  ImageToImageRGBTop->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  ImageToImageRGBTop->SetColormap( ColormapTop );

  ColormapBottom = Colormap::New();
  ColormapBottom->SetStandardColormap( PALETTE_GRAY );
  ColormapBottom->SetDataRange( 0, 2000 );
  ImageToImageRGBBottom = ImageToImageRGB::New();
  ImageToImageRGBBottom->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  ImageToImageRGBBottom->SetColormap( ColormapBottom );

  ColormapAlpha = Colormap::New();
  ColormapAlpha->SetStandardColormap( PALETTE_GRAY );
  ColormapAlpha->SetDataRange( 0, 2000 );
  ImageToImageRGBAlpha = ImageToImageRGB::New();
  ImageToImageRGBAlpha->SetAlphaMode( ImageToImageRGB::AlphaModeRamp );
  ImageToImageRGBAlpha->SetColormap( ColormapAlpha );

  FusionFilter = FusionAlpha::New();
  FusionFilter->SetInput( 0, ImageToImageRGBTop->GetOutput() );
  FusionFilter->SetInput( 1, ImageToImageRGBBottom->GetOutput() );
  FusionFilter->SetTransparencyImage( ImageToImageRGBAlpha->GetOutput() );
  FusionFilter->SetTopImageIndex( 0 );
  FusionFilter->SetOutputHasAlpha( true );

  View->slotConnectImage( FusionFilter->GetOutput() );

  SliderFrom = new QtSliderEntry( this );
  SliderFrom->slotSetTitle( "Transparent From" );
  QObject::connect( SliderFrom, SIGNAL( valueChanged( double ) ), this, SLOT( slotSliderFromChanged( double ) ) );
  ControlsLayout->addWidget( SliderFrom );
  SliderTo = new QtSliderEntry( this );
  SliderTo->slotSetTitle( "Transparent To" );
  QObject::connect( SliderTo, SIGNAL( valueChanged( double ) ), this, SLOT( slotSliderToChanged( double ) ) );
  ControlsLayout->addWidget( SliderTo );

  QObject::connect( this, SIGNAL( updateViewer() ), SLOT( slotUpdateSlice() ) );

  QObject::connect( this, SIGNAL( signalUpdate() ), SLOT( slotUpdateColormaps() ) );

  this->UpdateStudySelection();
  this->show();
}

QtFusionAlpha::~QtFusionAlpha()
{
  ImageToImageRGBTop->SetInput( NULL );
  ImageToImageRGBTop->Delete();
  ColormapTop->Delete();

  ImageToImageRGBBottom->SetInput( NULL );
  ImageToImageRGBBottom->Delete();
  ColormapBottom->Delete();

  FusionFilter->Delete();
}

void
QtFusionAlpha::UpdateStudySelection()
{
  // Get new list of study names.
  const QStringList *studyNamesList = FusionSlicer->GetTargetStudyNamesList();

  StudyNamesBoxTop->slotUpdateStudySelection( studyNamesList );
  this->slotSwitchStudyTop( StudyNamesBoxTop->GetCurrentName() );
  StudyNamesBoxBottom->slotUpdateStudySelection( studyNamesList );
  this->slotSwitchStudyBottom( StudyNamesBoxBottom->GetCurrentName() );
  StudyNamesBoxAlpha->slotUpdateStudySelection( studyNamesList );
  this->slotSwitchStudyAlpha( StudyNamesBoxAlpha->GetCurrentName() );
}

void
QtFusionAlpha::slotSliderFromChanged( double value )
{
  ImageToImageRGBAlpha->SetAlphaRampFrom( value );
  ImageToImageRGBAlpha->Update();
  View->slotRender();
}

void
QtFusionAlpha::slotSliderToChanged( double value )
{
  ImageToImageRGBAlpha->SetAlphaRampTo( value );
  ImageToImageRGBAlpha->Update();
  View->slotRender();
}

void
QtFusionAlpha::slotSwitchStudyTop( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.toLatin1() );

  if ( ! study.IsNull() ) 
    {
    StudyTop = study;
    ColormapTop->SetFromStudy( study );
    ImageToImageRGBTop->SetInput( FusionSlicer->GetOutput( study ) );
    
    if ( this->TopIsAlpha )
      {
      this->slotSwitchStudyAlpha( studyName );
      }
    else
      {
      View->slotRender();
      }
    }
}

void
QtFusionAlpha::slotSwitchStudyBottom( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.toLatin1() );

  if ( ! study.IsNull() ) 
    {
    StudyBottom = study;
    ColormapBottom->SetFromStudy( study );
    ImageToImageRGBBottom->SetInput( FusionSlicer->GetOutput( study ) );
    View->slotRender();
    }
}

void
QtFusionAlpha::slotSwitchStudyAlpha( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.toLatin1() );

  if ( ! study.IsNull() ) 
    {
    StudyAlpha = study;
    ColormapAlpha->SetFromStudy( study );
    
    float min = study->GetMinimumValue();
    float max = study->GetMaximumValue();
    const float range = max - min;
    const int precision = std::max( 0, -static_cast<int>( log( range / 256 ) / log( 10.0 ) ) );

    min -= range;
    max += range;

    SliderFrom->slotSetPrecision( precision );
    SliderFrom->slotSetRange( min, max );
    SliderTo->slotSetPrecision( precision );
    SliderTo->slotSetRange( min, max );

    ImageToImageRGBAlpha->SetInput( FusionSlicer->GetOutput( study ) );
    View->slotRender();
    }
}

void
QtFusionAlpha::slotUpdateSlice()
{
  View->GetRenderImage()->SetZoomFactorPercent( ZoomFactorPercent );
  View->slotRender();
}

void
QtFusionAlpha::slotUpdateColormap( Study::SmartPtr& study )
{
  bool update = false;

  if ( study == StudyTop ) 
    {
    ColormapTop->SetFromStudy( study );
    ColormapTop->Update();
    update = true;
    }
  if ( study == StudyBottom ) 
    {
    ColormapBottom->SetFromStudy( study );
    ColormapBottom->Update();
    update = true;
    }
  if ( study == StudyAlpha ) 
    {
    // nothing to do.
    }
  if ( update ) View->slotRender();
}

void
QtFusionAlpha::slotUpdateColormaps()
{
  ColormapTop->SetFromStudy( StudyTop );
  ColormapBottom->SetFromStudy( StudyBottom );
  ColormapBottom->SetFromStudy( StudyAlpha );
  View->slotRender();
}

void 
QtFusionAlpha::Export
( const QString& path, const QString& format, const QStringList* )
{
  QString filename;
  filename.sprintf( path.toLatin1(), "alp" );
  if ( format == QString::null || format == "PPM" )
    View->GetRenderImage()->WritePPM( filename.toLatin1() );
  else 
    {
    QPixmap pixmap = View->GetRenderImage()->GetPixmap();
    if ( ! pixmap.isNull() ) 
      {
      if ( !pixmap.save( filename, format.toLatin1() ) )
	QMessageBox::warning( this, "Save failed", "Error saving file" );
      }
    }
}

} // namespace cmtk
