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

#include "cmtkQtFusionEdge.h"

#include <qlayout.h>
#include <qmessagebox.h>
#include <QPixmap>

namespace
cmtk
{

QtFusionEdge::QtFusionEdge( QtSimpleFusionApp *const fusionApp, QWidget *const parent, Qt::WFlags flags )
  : QtFusionWindowTemplate( fusionApp, parent, flags ),
    StudyEdge( NULL ), StudyBG( NULL ),
    StudyEdgeImage( new Study )
{
  WindowLayout->setDirection( QBoxLayout::LeftToRight );

  this->setWindowTitle( "Edge Blending" );

  View = new QtScrollRenderView( this );
  ViewLayout->addWidget( View );

  this->m_FusionEdge = FusionEdge::New();
  View->slotConnectImage( this->m_FusionEdge->GetOutput() );

  StudyNamesBoxEdge = new QtStudyNamesBox( this );
  StudyNamesBoxEdge->slotSetLabel( "Edge Study:" );
  QObject::connect( StudyNamesBoxEdge, SIGNAL( signalActivated(const QString&)), this, SLOT( slotSwitchStudyEdge( const QString& ) ) );
  ControlsLayout->addWidget( StudyNamesBoxEdge );

  StudyNamesBoxBG = new QtStudyNamesBox( this );
  StudyNamesBoxBG->slotSetLabel( "Background Study:" );
  QObject::connect( StudyNamesBoxBG, SIGNAL( signalActivated(const QString&)), this, SLOT( slotSwitchStudyBG( const QString& ) ) );
  ControlsLayout->addWidget( StudyNamesBoxBG );

  SliderFrom = new QSlider( Qt::Horizontal, this );
  QObject::connect( SliderFrom, SIGNAL( valueChanged( int ) ), this, SLOT( slotSliderFromChanged( int ) ) );
  ControlsLayout->addWidget( SliderFrom );
  SliderTo = new QSlider( Qt::Horizontal, this );
  QObject::connect( SliderTo, SIGNAL( valueChanged( int ) ), this, SLOT( slotSliderToChanged( int ) ) );
  ControlsLayout->addWidget( SliderTo );

  WindowLevelBox = new QtWindowLevelControls( this );
  WindowLevelBox->slotSetStudy( StudyEdgeImage );
  QObject::connect( WindowLevelBox, SIGNAL( colormap( Study::SmartPtr& ) ), this, SLOT( slotUpdateColormap( Study::SmartPtr& ) ) );
  ControlsLayout->addWidget( WindowLevelBox );

  OperatorBox = new QComboBox( this );
  //  OperatorBox->setLabel( "Edge Operator" );
  OperatorBox->addItem( "Laplace" );
  OperatorBox->addItem( "Sobel" );
  ControlsLayout->addWidget( OperatorBox );

  QObject::connect( OperatorBox, SIGNAL( activated( int ) ), this, SLOT( slotParametersChanged() ) );

  SmoothCheckBox = new QCheckBox( this );
  SmoothCheckBox->setText( "Gaussian Smoothing" );
  ControlsLayout->addWidget( SmoothCheckBox );

  QObject::connect( SmoothCheckBox, SIGNAL( toggled( bool ) ), this, SLOT( slotParametersChanged() ) );
    
  GaussianWidthSlider = new QtSliderEntry( this );
  GaussianWidthSlider->slotSetTitle( "Gaussian Kernel Width" );
  GaussianWidthSlider->slotSetPrecision( 1 );
  GaussianWidthSlider->slotSetRange( 0.1, 4 );
  GaussianWidthSlider->slotSetValue( 0.5 );
  ControlsLayout->addWidget( GaussianWidthSlider );  

  QObject::connect( GaussianWidthSlider, SIGNAL( valueChanged( double ) ), this, SLOT( slotParametersChanged() ) );

  QObject::connect( this, SIGNAL( updateViewer() ), SLOT( slotUpdateSlice() ) );

  QObject::connect( this, SIGNAL( signalUpdate() ), SLOT( slotUpdateColormaps() ) );

  this->UpdateStudySelection();
  this->show();
}

QtFusionEdge::~QtFusionEdge()
{
  this->m_FusionEdge->Delete();
}

void
QtFusionEdge::UpdateStudySelection()
{
  // Get new list of study names.
  const QStringList *studyNamesList = FusionSlicer->GetTargetStudyNamesList();

  StudyNamesBoxEdge->slotUpdateStudySelection( studyNamesList );
  this->slotSwitchStudyEdge( StudyNamesBoxEdge->GetCurrentName() );
  StudyNamesBoxBG->slotUpdateStudySelection( studyNamesList );
  this->slotSwitchStudyBG( StudyNamesBoxBG->GetCurrentName() );
}

void
QtFusionEdge::slotSliderFromChanged( int value )
{
  this->m_FusionEdge->LookupEdge->SetAlphaRampFrom( value );
  View->slotRender();
}

void
QtFusionEdge::slotSliderToChanged( int value )
{
  this->m_FusionEdge->LookupEdge->SetAlphaRampTo( value );
  View->slotRender();
}

void
QtFusionEdge::slotSetColormap( int value )
{
  EdgeStandardColormap = value;
  StudyEdgeImage->SetStandardColormap( EdgeStandardColormap );
  this->slotUpdateColormap( StudyEdgeImage );
}

void
QtFusionEdge::slotParametersChanged()
{
  this->m_FusionEdge->EdgeOperator->SetOperator( OperatorBox->currentIndex() );
  this->m_FusionEdge->EdgeOperator->SetSmoothBeforeEdge( SmoothCheckBox->checkState() == Qt::Checked );
  this->m_FusionEdge->EdgeOperator->SetGaussianWidth( GaussianWidthSlider->GetValue() );
  View->slotRender();
}

void
QtFusionEdge::slotSwitchStudyEdge( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.toLatin1() );
  this->m_FusionEdge->SetInput( 0, FusionSlicer->GetOutput( study ) );
  SliderFrom->setRange( static_cast<int>( study->GetMinimumValue() ), static_cast<int>( study->GetMaximumValue() ) );
  SliderTo->setRange( static_cast<int>( study->GetMinimumValue() ), static_cast<int>( study->GetMaximumValue() ) );

  if ( study != StudyEdgeImage ) 
    {
    StudyEdgeImage->CopyColormap( study );
    StudyEdgeImage->SetStandardColormap( EdgeStandardColormap );
    WindowLevelBox->slotSetStudy( StudyEdgeImage );
    this->slotUpdateColormap( StudyEdgeImage );
    }
  
  View->slotRender();
}

void
QtFusionEdge::slotSwitchStudyBG( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.toLatin1() );
  this->m_FusionEdge->SetInput( 1, FusionSlicer->GetOutput( study ) );
  View->slotRender();

  if ( study != StudyBG ) 
    {
    StudyBG = study;
    this->slotUpdateColormap( StudyBG );
    }
}

void
QtFusionEdge::slotUpdateSlice()
{
  View->GetRenderImage()->SetZoomFactorPercent( ZoomFactorPercent );
  View->slotRender();
}

void
QtFusionEdge::slotUpdateColormap( Study::SmartPtr& study )
{
  bool update = false;

  if ( study == StudyBG ) 
    {
    this->m_FusionEdge->ColormapImage->SetFromStudy( study );
    update = true;
  }
  if ( study == StudyEdgeImage ) 
    {
    this->m_FusionEdge->ColormapEdge->SetFromStudy( study );
    update = true;
    }
  
  if ( update ) 
    View->slotRender();
}

void
QtFusionEdge::slotUpdateColormaps()
{
  this->m_FusionEdge->ColormapImage->SetFromStudy( StudyBG );
  View->slotRender();
}

void 
QtFusionEdge::Export
( const QString& path, const QString& format, const QStringList* )
{
  QString filename;
  filename.sprintf( path.toLatin1(), "edg" );

  QString fmt = format;
  if ( fmt == QString::null )
    fmt = path.section( '.', -1 ).toUpper();

  QPixmap pixmap = View->GetRenderImage()->GetPixmap();
  if ( ! pixmap.isNull() ) 
    {
    if ( !pixmap.save( filename, fmt.toLatin1() ) )
      QMessageBox::warning( this, "Save failed", "Error saving file" );
    }
}

} // namespace cmtk
