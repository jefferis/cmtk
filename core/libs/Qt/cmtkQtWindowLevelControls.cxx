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

#include <cmtkQtWindowLevelControls.h>

#include <cmtkColormap.h>

#include <qcombobox.h>
#include <Q3VBoxLayout>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtWindowLevelControls::QtWindowLevelControls
( QWidget *const parent, const char* name )
  : QWidget( parent, name ),
    Study( NULL )
{
  Layout = new Q3VBoxLayout( this );
  Layout->setContentsMargins( 5, 5, 5, 5 );

  QComboBox* colormapBox = new QComboBox( this );
  Layout->addWidget( colormapBox );
  
  for ( unsigned int colormapIndex = 0; Colormap::StandardColormaps[colormapIndex]; ++colormapIndex ) 
    {
    colormapBox->insertItem( Colormap::StandardColormaps[colormapIndex] );
    }
  
  QObject::connect( colormapBox, SIGNAL( activated( int ) ), this, SLOT( slotSelectColormap( int ) ) );

  BlackWindowSlider = new QtSliderEntry( this );
  QObject::connect( BlackWindowSlider, SIGNAL( valueChanged( double ) ), this, SLOT( slotControlsChanged() ) );
  BlackWindowSlider->slotSetTitle( "Black" );
  BlackWindowSlider->slotSetMinMaxLabels( QString::null, QString::null );
  Layout->addWidget( BlackWindowSlider );

  WhiteLevelSlider = new QtSliderEntry( this );
  QObject::connect( WhiteLevelSlider, SIGNAL( valueChanged( double ) ), this, SLOT( slotControlsChanged() ) );
  WhiteLevelSlider->slotSetTitle( "White" );
  WhiteLevelSlider->slotSetMinMaxLabels( QString::null, QString::null );
  Layout->addWidget( WhiteLevelSlider );

  WindowLevelCheckBox = new QCheckBox( "Window/Level", this );
  QObject::connect( WindowLevelCheckBox, SIGNAL( stateChanged( int ) ), this, SLOT( slotSwitchModeWL( int ) ) );
  Layout->addWidget( WindowLevelCheckBox );

  GammaSlider = new QtSliderEntry( this );
  GammaSlider->slotSetPrecision( 1 );
  GammaSlider->slotSetRange( 0.1, 10 );
  GammaSlider->slotSetValue( 1 );
  GammaSlider->slotSetTitle( "Gamma Value" );
  GammaSlider->slotSetMinMaxLabels( QString::null, QString::null );
  QObject::connect( GammaSlider, SIGNAL( valueChanged( double ) ), this, SLOT( slotControlsChanged() ) );
  Layout->addWidget( GammaSlider );

  Layout->addItem( new QSpacerItem( 0, 0, QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding ) );
}

void
QtWindowLevelControls::slotSetStudy( Study::SmartPtr& study )
{
  Study = study;
  RangeFrom = Study->GetMinimumValue();
  RangeTo = Study->GetMaximumValue();

  RangeWidth = RangeTo - RangeFrom;

  const Volume* volume = Study->GetVolume();
  if ( volume ) 
    {
    const TypedArray* data = volume->GetData();
    if ( data ) 
      {
      Types::DataItem mean, variance;
      data->GetStatistics( mean, variance );
      RangeWidth = sqrt( variance );
      }
    }
  
  this->slotSwitchModeWL( WindowLevelCheckBox->isChecked() );
}

void 
QtWindowLevelControls::slotSwitchModeWL( int modeWindowLevel )
{
  if ( Study.IsNull() ) return;
  
  const float black = Study->GetBlack();
  const float white = Study->GetWhite();
  
  unsigned int precision = 0;
  
  if ( RangeWidth > 0 ) 
    {
    precision = static_cast<unsigned int>( std::max( 0.0, (log( 1.0 / 256 ) + log( RangeWidth )) / log(0.1) ) );
    }
  
  WhiteLevelSlider->slotSetPrecision( precision );
  BlackWindowSlider->slotSetPrecision( precision );

  if ( modeWindowLevel ) 
    {
    BlackWindowSlider->slotSetRange( 0, RangeTo - RangeFrom );
    BlackWindowSlider->slotSetValue( white - black );
    BlackWindowSlider->slotSetTitle( "Window" );
    
    WhiteLevelSlider->slotSetRange( RangeFrom, RangeTo );
    WhiteLevelSlider->slotSetValue( (white + black) / 2 );
    WhiteLevelSlider->slotSetTitle( "Level" );
    } 
  else
    {
    BlackWindowSlider->slotSetRange( RangeFrom, RangeTo );
    BlackWindowSlider->slotSetValue( black );
    BlackWindowSlider->slotSetTitle( "Black" );
    
    WhiteLevelSlider->slotSetRange( RangeFrom, RangeTo );
    WhiteLevelSlider->slotSetValue( white );
    WhiteLevelSlider->slotSetTitle( "White" );
    }
}

void
QtWindowLevelControls::slotControlsChanged()
{
  if ( Study.IsNull() ) return;

  float black, white;
  if ( WindowLevelCheckBox->isChecked() ) 
    {
    black = WhiteLevelSlider->GetValue() - BlackWindowSlider->GetValue() / 2;
    white = WhiteLevelSlider->GetValue() + BlackWindowSlider->GetValue() / 2;
    } 
  else
    {
    black = BlackWindowSlider->GetValue();
    white = WhiteLevelSlider->GetValue();
    }
  
  float gamma = GammaSlider->GetValue();

  Study->SetBlack( black );
  Study->SetWhite( white );
  Study->SetGamma( gamma );

  emit colormap( Study );
}

void
QtWindowLevelControls::slotSelectColormap( int colormapIndex )
{
  if ( !Study.IsNull() )
    {
    Study->SetStandardColormap( colormapIndex );
    emit colormap( Study );
    }
}

} // namespace cmtk
