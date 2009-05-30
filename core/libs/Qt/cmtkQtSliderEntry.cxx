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

#include <cmtkQtSliderEntry.h>

#include <stdlib.h>
#include <math.h>

#include <qfont.h>

#include <QGridLayout>
#include <QLabel>

#include <cmtkMathUtil.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtSliderEntry::QtSliderEntry( QWidget* parent, const char *name )
  : QWidget( parent, name )
{
  QFont font = this->font();
  font.setPointSize( 2 * font.pointSize() / 4 );
  
  Layout = new QGridLayout( this, 3, 4 );
  Layout->setColStretch( 0, 1 );
  Layout->setColStretch( 1, 1 );
  Layout->setColStretch( 2, 0 );
  Layout->setColStretch( 3, 0 );

  Slider = new QSlider( Qt::Horizontal, this );
  QObject::connect( Slider, SIGNAL( valueChanged( int ) ), this, SLOT( slotSliderValueChanged( int ) ) );
  Layout->addMultiCellWidget( Slider, 1, 1, 0, 1 );

  Edit = new QLineEdit( this );
  Edit->setFixedWidth( 100 );
  Validator = new QDoubleValidator( Edit );
  Edit->setValidator( Validator );
  QObject::connect( Edit, SIGNAL( returnPressed() ), this, SLOT( slotEditReturnPressed() ) );
  Layout->addWidget( Edit, 1, 3 );

  TitleLabel = new QLabel( this );
  TitleLabel->hide();
  MinLabel = new QLabel( this );
  MinLabel->setFont( font );
  MinLabel->hide();
  MaxLabel = new QLabel( this );
  MaxLabel->setFont( font );
  MaxLabel->setAlignment( Qt::AlignRight );
  MaxLabel->hide();

  Precision = 0;
  PrecisionFactor = 1;
}

double 
QtSliderEntry::GetValue() const
{
  return Edit->text().toDouble();
}

double
QtSliderEntry::GetMinValue() const
{
  return 1.0 * Slider->minValue() / PrecisionFactor;
}

double
QtSliderEntry::GetMaxValue() const
{
  return 1.0 * Slider->minValue() / PrecisionFactor;
}

void
QtSliderEntry::slotSetRange( double rangeFrom, double rangeTo )
{
  const double rangeWidth = rangeTo - rangeFrom;
  
  if ( rangeWidth > 0 ) 
    {
    const int autoPrecision = static_cast<int>( (log( 0.001 ) + log( rangeWidth )) / log( 0.1 ) );
    this->slotSetPrecision( std::max<int>( this->Precision, autoPrecision ) );
    }
  
  Slider->setRange( static_cast<int>( rangeFrom * PrecisionFactor ), static_cast<int>( rangeTo * PrecisionFactor ) );
  
  Validator->setRange( rangeFrom - 10 * rangeWidth, rangeTo + 10 * rangeWidth, Precision );
  
  MinLabel->setNum( rangeFrom );
  MaxLabel->setNum( rangeTo );
}

void
QtSliderEntry::slotSetPrecision( int precision )
{
  Precision = precision;
  PrecisionFactor = static_cast<uint>( pow( 10.0, precision ) );
  Validator->setDecimals( Precision );
}

void
QtSliderEntry::slotSetTitle( const QString& title )
{
  TitleLabel->setText( title );
  Layout->addMultiCellWidget( TitleLabel, 0, 0, 0, 2 );
  TitleLabel->show();
}

void
QtSliderEntry::slotSetMinMaxLabels( const QString& minLabel, const QString& maxLabel )
{
  if ( minLabel != QString::null ) 
    {
    MinLabel->setText( minLabel );
    } 
  else
    {
    MinLabel->setNum( Validator->bottom() );
    }
  Layout->addWidget( MinLabel, 2, 0 );
  MinLabel->show();
  
  if ( maxLabel != QString::null ) 
    {
    MaxLabel->setText( maxLabel );
    } 
  else
    {
    MaxLabel->setNum( Validator->top() );
    }
  Layout->addWidget( MaxLabel, 2, 1 );
  MaxLabel->show();
}

void 
QtSliderEntry::slotEditReturnPressed()
{
  double value = atof( Edit->text() );
  int valueSlider = static_cast<int>( value * PrecisionFactor );

  if ( valueSlider < Slider->minValue() )
    this->slotSetRange( value, Slider->maxValue() / PrecisionFactor );
  if ( valueSlider > Slider->maxValue() )
    this->slotSetRange( Slider->minValue() / PrecisionFactor, value );

  Slider->setValue( valueSlider );
  emit valueChanged( value );
}

void 
QtSliderEntry::slotSliderValueChanged( int value )
{
  double realValue = 1.0 * value / PrecisionFactor;
  QString valueString;
  Edit->setText( valueString.setNum( realValue, 'f', Precision ) );
  emit valueChanged( realValue );
}

void 
QtSliderEntry::slotCenter()
{
  Slider->setValue( (Slider->minValue() + Slider->maxValue()) / 2 );
  // valueChanged signal should be emitted indirectly by slotSliderValueChanged
}

void 
QtSliderEntry::slotSetValue( const double value )
{
  QString valueString;
  Edit->setText( valueString.setNum( value, 'f', Precision ) );
  int valueSlider = static_cast<int>( value * PrecisionFactor );

  if ( valueSlider < Slider->minValue() )
    this->slotSetRange( value, Slider->maxValue() / PrecisionFactor );
  if ( valueSlider > Slider->maxValue() )
    this->slotSetRange( Slider->minValue() / PrecisionFactor, value );

  Slider->setValue( valueSlider );
  emit valueChanged( value );
}

} // namespace cmtk
