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

#ifndef __cmtkQtSliderEntry_h_included_
#define __cmtkQtSliderEntry_h_included_

#include <cmtkconfig.h>

#include <qwidget.h>
#include <qslider.h>
#include <qlabel.h>
#include <qlineedit.h>
#include <qvalidator.h>
#include <qlayout.h>

#include <QGridLayout>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/** Widget that combines a slider with a numerical entry field and labels.
 */
class QtSliderEntry :
  /// We use a vertical group box as the base class.
  public QWidget
{
  Q_OBJECT // we use slots and signals.

public:
  /// Constructor.
  QtSliderEntry( QWidget* parent, const char* name = NULL );

  /// Get value.
  double GetValue() const;

  /// Get minimum value.
  double GetMinValue() const;

  /// Get maximum value.
  double GetMaxValue() const;

signals:
  /// Emitted when value changes.
  void valueChanged( double value );

public slots:
  /// Set title label.
  void slotSetTitle( const QString& title );

  /// Set min/max labels.
  void slotSetMinMaxLabels( const QString& minLabel, const QString& maxLabel );

  /// Set value range.
  void slotSetRange( double rangeFrom, double rangeTo );

  /// Set number of digits.
  void slotSetPrecision( int precision );

  /// Set value.
  void slotSetValue( const double value );

  /// Set to center position.
  void slotCenter();

private slots:
  /// Called when "Return" is pressed in the line edit field.
  void slotEditReturnPressed();

  /// Called when line edit value changes.
  void slotSliderValueChanged( int value );

private:
  /// Number of decimal digits.
  uint Precision;

  /// Factor to convert between integer and true float representation.
  uint PrecisionFactor;
  
  /// Layout for children.
  QGridLayout* Layout;

  /// The slider object.
  QSlider* Slider;
  
  /// The entry field.
  QLineEdit* Edit;

  /// The entry validator object.
  QDoubleValidator* Validator;

  /// The label for the widget title.
  QLabel* TitleLabel;

  /// The label for the slider's minimum value.
  QLabel* MinLabel;

  /// The label for the slider's maximum value.
  QLabel* MaxLabel;

};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtSliderEntry_h_included_
