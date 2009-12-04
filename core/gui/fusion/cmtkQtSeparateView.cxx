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

#include <cmtkQtSeparateView.h>

#include <qlayout.h>
#include <QVBoxLayout>

namespace
cmtk
{

QtSeparateView::QtSeparateView( QtSimpleFusionApp* fusionApp, QWidget *const parent, Qt::WFlags flags )
  : QtFusionWindowTemplate( fusionApp, parent, "SeparateViewWindow", flags )
{
  this->setCaption( "Separate Display" );

  // transform the default vertical view layout into a horizontal layout
  ViewLayout->setDirection( QBoxLayout::LeftToRight );

  this->Left.Construct( this, ViewLayout, "Left:" );
  QObject::connect( this->Left.StudyNamesBox, SIGNAL( signalActivated(const QString&) ), this, SLOT( slotSwitchStudyL( const QString& ) ) );

  this->Right.Construct( this, ViewLayout, "Right:" );
  QObject::connect( this->Right.StudyNamesBox, SIGNAL( signalActivated(const QString&) ), this, SLOT( slotSwitchStudyR( const QString& ) ) );

  QObject::connect( this, SIGNAL( signalUpdate() ), SLOT( slotUpdateColormaps() ) );
  QObject::connect( this, SIGNAL( updateViewer() ), SLOT( slotUpdateSlice() ) );

  this->UpdateStudySelection();
  this->show();
}

void
QtSeparateView::UI::Construct
( QWidget *const parent, QLayout *const inLayout, const QString& label )
{
  QVBoxLayout* layout = new QVBoxLayout( inLayout );
  StudyNamesBox = new QtStudyNamesBox( parent, "StudyNamesBox" );
  StudyNamesBox->slotSetLabel( label );
  layout->addWidget( StudyNamesBox );
  View = new QtScrollRenderView( parent );
  layout->addWidget( View );

  this->m_Colormap = Colormap::New();
  this->m_Colormap->SetStandardColormap( PALETTE_GRAY );
  this->m_Colormap->SetDataRange( 0, 2000 );
  this->m_ImageToImageRGB = ImageToImageRGB::New();
  this->m_ImageToImageRGB->SetAlphaMode( ImageToImageRGB::AlphaModeConst );
  this->m_ImageToImageRGB->SetColormap( this->m_Colormap );
  View->slotConnectImage( this->m_ImageToImageRGB->GetOutput() );
}

void
QtSeparateView::UpdateStudySelection()
{
  // Get new list of study names.
  const QStringList *studyNamesList = FusionSlicer->GetTargetStudyNamesList();

  this->Left.StudyNamesBox->slotUpdateStudySelection( studyNamesList );
  this->slotSwitchStudyL( this->Left.StudyNamesBox->GetCurrentName() );
  this->Right.StudyNamesBox->slotUpdateStudySelection( studyNamesList );
  this->slotSwitchStudyR( this->Right.StudyNamesBox->GetCurrentName() );
}

void
QtSeparateView::slotUpdateColormaps()
{
  if ( this->Left.m_Study ) 
    {
    this->Left.m_Colormap->SetFromStudy( this->Left.m_Study );
    this->Left.View->slotRender();
    }
  if ( this->Right.m_Study ) 
    {
    this->Right.m_Colormap->SetFromStudy( this->Right.m_Study );
    this->Right.View->slotRender();
    }
}

void
QtSeparateView::slotSwitchStudyL( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.latin1() );
  if ( study ) 
    {
    this->Left.m_Study = study;
    this->Left.m_Colormap->SetFromStudy( this->Left.m_Study );
    this->Left.m_ImageToImageRGB->SetInput( FusionSlicer->GetOutput( this->Left.m_Study ) );
    this->Left.View->slotRender();
    }
}

void
QtSeparateView::slotSwitchStudyR( const QString& studyName )
{
  Study::SmartPtr study = FusionApp->m_StudyList->FindStudyName( studyName.latin1() );
  if ( study ) 
    {
    this->Right.m_Study = study;
    this->Right.m_Colormap->SetFromStudy( this->Right.m_Study );
    this->Right.m_ImageToImageRGB->SetInput( FusionSlicer->GetOutput( study ) );
    this->Right.View->slotRender();
    }
}

void
QtSeparateView::slotUpdateSlice()
{
  this->Left.View->GetRenderImage()->SetZoomFactorPercent( ZoomFactorPercent );
  this->Left.View->slotRender();
  this->Right.View->GetRenderImage()->SetZoomFactorPercent( ZoomFactorPercent );
  this->Right.View->slotRender();
}

void
QtSeparateView::slotUpdateColormap( Study::SmartPtr& study )
{
  if ( study == this->Left.m_Study ) 
    this->Left.m_Colormap->SetFromStudy( study );
  this->Left.View->slotRender();

  if ( study == this->Right.m_Study )
    this->Right.m_Colormap->SetFromStudy( study );
  this->Right.View->slotRender();
}

} // namespace cmtk
