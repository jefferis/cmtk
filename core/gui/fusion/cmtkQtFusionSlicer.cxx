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

#include "cmtkQtFusionSlicer.h"

#include "cmtkQtFusionGlobal.h"
#include "cmtkQtSimpleFusionApp.h"

#include "Pipeline/cmtkRendererCollection.h"

#include <qlayout.h>
#include <qgroupbox.h>
#include <qbuttongroup.h>
#include <qradiobutton.h>
#include <qcheckbox.h>
#include <qpushbutton.h>

#include <QVBoxLayout>

namespace
cmtk
{

QtFusionSlicer::QtFusionSlicer( QtSimpleFusionApp *const fusionApp ) :
  QWidget( NULL, 0 ),
  FusionApp( fusionApp )
{
  this->setWindowIcon( QtFusionGlobal::WindowIcon() );
  this->setWindowTitle( "Planar Slicer" );
  this->hide();
  this->setSizePolicy( QSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed ) );

  this->m_FusionSlicers = new FusionSlicers;
  StudyNamesList = new QStringList;
  TargetStudyNamesList = new QStringList;

  QVBoxLayout* layout = new QVBoxLayout( this );

  ReferenceBox = new QtStudyNamesBox( this );
  ReferenceBox->slotSetLabel( "Reference Study:" );
  layout->addWidget( ReferenceBox, 1 );
  QObject::connect( ReferenceBox, SIGNAL( signalActivated( const QString& ) ), this, SLOT( slotSetReferenceStudy( const QString& ) ) );

  QGroupBox* orientButtonBox = new QGroupBox( "Orientation", this );

  QButtonGroup* orientButtonGroup = new QButtonGroup( orientButtonBox );
  orientButtonGroup->setExclusive( true );
  
  QVBoxLayout *orientLayout = new QVBoxLayout;
  orientButtonBox->setLayout( orientLayout );

  QRadioButton* axialButton = new QRadioButton( "Axial", orientButtonBox );
  orientLayout->addWidget( axialButton );
  orientButtonGroup->addButton( axialButton, AXIS_Z );
  
  QRadioButton* sagittalButton = new QRadioButton( "Sagittal", orientButtonBox );
  orientLayout->addWidget( sagittalButton );
  orientButtonGroup->addButton( sagittalButton, AXIS_X );

  QRadioButton* coronalButton = new QRadioButton( "Coronal", orientButtonBox );
  orientLayout->addWidget( coronalButton );
  orientButtonGroup->addButton( coronalButton, AXIS_Y );

  axialButton->setChecked( true );

  QObject::connect( orientButtonGroup, SIGNAL( buttonClicked( int ) ), this, SLOT( slotSetOrientation( int ) ) );
  layout->addWidget( orientButtonBox, 1 );

  QGroupBox* sliceGroupBox = new QGroupBox( this );
  sliceGroupBox->setTitle( "Slice Plane Parameters" );
  layout->addWidget( sliceGroupBox );

  SliceSlider = new QtSliderEntry( sliceGroupBox );
  SliceSlider->slotSetPrecision( 2 );
  QObject::connect( SliceSlider, SIGNAL( valueChanged( double ) ), this, SLOT( slotSetSlicePosition( double ) ) );

  QPushButton* centerButton = new QPushButton( "Center", sliceGroupBox );
  QObject::connect( centerButton, SIGNAL( clicked() ), SliceSlider, SLOT( slotCenter() ) );

  QVBoxLayout *sliceLayout = new QVBoxLayout;
  sliceLayout->addWidget( SliceSlider );
  sliceLayout->addWidget( centerButton );
  sliceGroupBox->setLayout( sliceLayout );
  
  QGroupBox* interButtonBox = new QGroupBox( "Interpolation", this );

  QButtonGroup* interButtonGroup = new QButtonGroup( interButtonBox );
  interButtonGroup->setExclusive( true );
  
  QVBoxLayout *interLayout = new QVBoxLayout;
  interButtonBox->setLayout( interLayout );

  QRadioButton* nnButton = new QRadioButton( "Nearest Neighbor", interButtonBox );
  interLayout->addWidget( nnButton );
  interButtonGroup->addButton( nnButton, cmtk::Interpolators::NEAREST_NEIGHBOR );
  
  QRadioButton* linearButton = new QRadioButton( "Linear", interButtonBox );
  interLayout->addWidget( linearButton );
  interButtonGroup->addButton( linearButton, cmtk::Interpolators::LINEAR );
  
  QRadioButton* cubicButton = new QRadioButton( "Cubic", interButtonBox );
  interLayout->addWidget( cubicButton );
  interButtonGroup->addButton( cubicButton, cmtk::Interpolators::CUBIC );

  QRadioButton* pvButton = new QRadioButton( "Partial Volume", interButtonBox );
  interLayout->addWidget( pvButton );
  interButtonGroup->addButton( pvButton, cmtk::Interpolators::PARTIALVOLUME );
  
  linearButton->setChecked( true );

  QObject::connect( interButtonGroup, SIGNAL( buttonClicked( int ) ), this, SLOT( slotSetInterpolation( int ) ) );
  layout->addWidget( interButtonBox );

  QCheckBox* warpBox = new QCheckBox( "Apply Warp", this );
  warpBox->setChecked( false );
  layout->addWidget( warpBox );
  QObject::connect( warpBox, SIGNAL( stateChanged( int ) ), this, SLOT( slotApplyWarp( int ) ) );  

  QPushButton* closeButton = new QPushButton( "Close", this );
  QObject::connect( closeButton, SIGNAL( clicked() ), this, SLOT( hide() ) );
  layout->addWidget( closeButton );

  QObject::connect( FusionApp, SIGNAL( signalStudyListChanged() ), this, SLOT( slotStudyListChanged() ) );
  QObject::connect( FusionApp, SIGNAL( signalReferenceStudyChanged() ), this, SLOT( slotReferenceStudyChanged() ) );

  this->slotSetOrientation( AXIS_Z );
  this->slotSetInterpolation( cmtk::Interpolators::LINEAR );

  this->slotStudyListChanged();
  this->slotReferenceStudyChanged();
}

QtFusionSlicer::~QtFusionSlicer()
{
  if ( this->m_FusionSlicers ) 
    this->m_FusionSlicers->Delete();
  if ( StudyNamesList ) 
    delete StudyNamesList;
  if ( TargetStudyNamesList ) 
    delete TargetStudyNamesList;
}

void
QtFusionSlicer::slotStudyListChanged()
{
  StudyNamesList->clear();

  if ( FusionApp->m_StudyList ) 
    {
    const unsigned int numberOfStudies = FusionApp->m_StudyList->size();
    for ( unsigned int idx = 0; idx < numberOfStudies; ++idx ) 
      {
      QString newStudyName( FusionApp->m_StudyList->GetStudy( idx )->GetName() );
      if ( ! StudyNamesList->contains( newStudyName ) )
	StudyNamesList->push_back( newStudyName );
      }
    this->m_FusionSlicers->SetStudyList( FusionApp->m_StudyList );
    }
  
  ReferenceBox->slotUpdateStudySelection( StudyNamesList );
}

void
QtFusionSlicer::slotSetReferenceStudy( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.toLatin1() );
  if ( study ) 
    {
    FusionApp->slotSetReferenceStudy( study );
    }
}

void
QtFusionSlicer::slotReferenceStudyChanged()
{
  if ( FusionApp->ReferenceStudy ) 
    {
    FusionApp->ReferenceStudy->ReadVolume( false /*reread*/, AnatomicalOrientation::ORIENTATION_STANDARD );
    this->m_FusionSlicers->SetReferenceStudy( FusionApp->ReferenceStudy );
    ReferenceBox->slotSetCurrentText( FusionApp->ReferenceStudy->GetName() );
    this->slotSetOrientation( AXIS_Z );

    FusionApp->m_StudyList->Print();

    TargetStudyNamesList->clear();
    TargetStudyNamesList->push_back( FusionApp->ReferenceStudy->GetName() );

    StudyToXform::const_iterator it = (*FusionApp->m_StudyList)[FusionApp->ReferenceStudy].begin();
    while ( it != (*FusionApp->m_StudyList)[FusionApp->ReferenceStudy].end() ) 
      {
      QString name( it->first->GetName() );
      if ( ! TargetStudyNamesList->contains( name ) ) 
	TargetStudyNamesList->push_back( name );      
      ++it;
      }
    }
}

void 
QtFusionSlicer::slotSetOrientation( int sliceNormal )
{
  this->SliceNormal = sliceNormal;
  switch ( SliceNormal ) 
    {
    case 0:
      SliceSlider->slotSetMinMaxLabels( "Left", "Right" );
      break;
    case 1:
      SliceSlider->slotSetMinMaxLabels( "Ventral", "Dorsal" );
      break;
    case 2:
      SliceSlider->slotSetMinMaxLabels( "Caudal", "Cranial" );
      break;
    }
  this->m_FusionSlicers->SetSliceNormal( SliceNormal );
  
  if ( ! FusionApp->ReferenceStudy ) return;
  const UniformVolume *volume = FusionApp->ReferenceStudy->GetVolume();
  if ( ! volume ) return;
  
  
  this->PositionTo = volume->Size[SliceNormal];
  this->PositionStep = volume->m_Delta[SliceNormal];
  
  SliceSlider->slotSetRange( 0, this->PositionTo );
  SliceSlider->slotCenter();
}

void 
QtFusionSlicer::slotSetSlicePosition( double slicePosition )
{
  if ( this->m_FusionSlicers ) 
    {
    this->m_FusionSlicers->SetReferenceSlicePosition( slicePosition );
    emit sliceChanged();
    }
}

void
QtFusionSlicer::slotSetInterpolation( int mode )
{
  if ( this->m_FusionSlicers ) 
    {
    this->m_FusionSlicers->SetInterpolationMode( (cmtk::Interpolators::InterpolationEnum) mode );
    emit sliceChanged();
    }
}

void
QtFusionSlicer::slotApplyWarp( int mode )
{
  if ( this->m_FusionSlicers ) 
    {
    this->m_FusionSlicers->SetApplyWarp( mode );
    emit sliceChanged();
    }
}

Image* 
QtFusionSlicer::GetOutput( Study::SmartPtr& study )
{
  return this->m_FusionSlicers->GetOutput( study );
}

} // namespace cmtk
