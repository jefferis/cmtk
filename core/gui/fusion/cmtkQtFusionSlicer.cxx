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

#include <cmtkQtFusionSlicer.h>

#include <cmtkQtFusionGlobal.h>

#include <qlayout.h>
#include <q3vgroupbox.h>
#include <qradiobutton.h>
#include <qcheckbox.h>
#include <qpushbutton.h>

#include <Q3VBoxLayout>
#include <Q3ButtonGroup>

#include <cmtkQtSimpleFusionApp.h>
#include <cmtkRendererCollection.h>

namespace
cmtk
{

QtFusionSlicer::QtFusionSlicer( QtSimpleFusionApp *const fusionApp ) :
  QWidget( NULL, "FusionSlicer", 0 ),
  FusionApp( fusionApp )
{
  this->setIcon( QtFusionGlobal::WindowIcon() );
  this->setCaption( "Planar Slicer" );
  this->hide();
  this->setSizePolicy( QSizePolicy( QSizePolicy::Fixed, QSizePolicy::Fixed ) );

  this->m_FusionSlicers = new FusionSlicers;
  StudyNamesList = new QStringList;
  TargetStudyNamesList = new QStringList;

  Q3VBoxLayout* layout = new Q3VBoxLayout( this );

  ReferenceBox = new QtStudyNamesBox( this, "ReferenceBox" );
  ReferenceBox->slotSetLabel( "Reference Study:" );
  layout->addWidget( ReferenceBox );
  QObject::connect( ReferenceBox, SIGNAL( signalActivated( const QString& ) ), this, SLOT( slotSetReferenceStudy( const QString& ) ) );

  Q3ButtonGroup* orientButtonGroup = new Q3VButtonGroup( "Orientation", this, "OrientButtonGroup" );
  orientButtonGroup->insert( new QRadioButton( "Axial", orientButtonGroup ), AXIS_Z );
  orientButtonGroup->insert( new QRadioButton( "Sagittal", orientButtonGroup ), AXIS_X );
  orientButtonGroup->insert( new QRadioButton( "Coronal", orientButtonGroup ), AXIS_Y );
  orientButtonGroup->setExclusive( true );
  orientButtonGroup->setButton( AXIS_Z );
  QObject::connect( orientButtonGroup, SIGNAL( clicked( int ) ), this, SLOT( slotSetOrientation( int ) ) );
  layout->addWidget( orientButtonGroup );

  Q3GroupBox* sliceGroupBox = new Q3VGroupBox( this, "SliceGroupBox" );
  sliceGroupBox->setTitle( "Slice Plane Parameters" );
  layout->addWidget( sliceGroupBox );

  SliceSlider = new QtSliderEntry( sliceGroupBox, "SliceSlider" );
  SliceSlider->slotSetPrecision( 2 );
  QObject::connect( SliceSlider, SIGNAL( valueChanged( double ) ), this, SLOT( slotSetSlicePosition( double ) ) );

  QPushButton* centerButton = new QPushButton( "Center", sliceGroupBox, "CenterButton" );
  QObject::connect( centerButton, SIGNAL( clicked() ), SliceSlider, SLOT( slotCenter() ) );
  
  Q3ButtonGroup* interButtonGroup = new Q3VButtonGroup( "Interpolation", this, "InterButtonGroup" );
  interButtonGroup->insert( new QRadioButton( "Nearest Neighbor", interButtonGroup ), cmtk::Interpolators::NEAREST_NEIGHBOR );
  interButtonGroup->insert( new QRadioButton( "Linear", interButtonGroup ), cmtk::Interpolators::LINEAR );
  interButtonGroup->insert( new QRadioButton( "Cubic", interButtonGroup ), cmtk::Interpolators::CUBIC );
  interButtonGroup->insert( new QRadioButton( "Partial Volume", interButtonGroup ), cmtk::Interpolators::PARTIALVOLUME );
  interButtonGroup->setExclusive( true );
  QObject::connect( interButtonGroup, SIGNAL( clicked( int ) ), this, SLOT( slotSetInterpolation( int ) ) );
  interButtonGroup->setButton( cmtk::Interpolators::LINEAR );
  layout->addWidget( interButtonGroup );

  QCheckBox* warpBox = new QCheckBox( "Apply Warp", this, "WarpCheckBox" );
  warpBox->setChecked( false );
  layout->addWidget( warpBox );
  QObject::connect( warpBox, SIGNAL( stateChanged( int ) ), this, SLOT( slotApplyWarp( int ) ) );
  

  QPushButton* closeButton = new QPushButton( "Close", this, "CloseButton" );
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
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.latin1() );
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
  this->PositionStep = volume->Delta[SliceNormal];
  
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
