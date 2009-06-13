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

#ifndef __cmtkQtWindowLevelControls_h_included_
#define __cmtkQtWindowLevelControls_h_included_

#include <cmtkconfig.h>

#include <q3groupbox.h>
#include <qcheckbox.h>
#include <qslider.h>
#include <qlayout.h>
#include <Q3VBoxLayout>

#include <cmtkQtSliderEntry.h>
#include <cmtkStudy.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/// Widget for a group box with Window/Level controls.
class QtWindowLevelControls :
  /// Inherit from Qt's group box.
  public QWidget
{
  Q_OBJECT // we're using signals and slots

signals:
  /// This signal is emitted when the controls change.
  void colormap( Study::SmartPtr& );

public slots:
  /// This signal tells the controls that the study object has changed.
  void slotSetStudy( Study::SmartPtr& study );

private slots:
  /// This slot is called when the Window/Level mode is changed.
  void slotSwitchModeWL( int );

  /// This slot is called by the UI widgets when their values change.
  void slotControlsChanged();

  /// This slot is called when the user picks a new colormap.
  void slotSelectColormap( int colormapIndex );

public:
  /// Constructor.
  QtWindowLevelControls( QWidget *const parent, const char* name = 0 );

private:
  /// The study object that we're working on.
  Study::SmartPtr m_Study;

  /// Layout of this widget.
  Q3VBoxLayout* Layout;

  /// The top slider in the UI.
  QtSliderEntry* BlackWindowSlider;

  /// The bottom slider in the UI.
  QtSliderEntry* WhiteLevelSlider;

  /// The gamma slider in the UI.
  QtSliderEntry* GammaSlider;

  /// The checkbox that switches between Window/Level and Black/White mode.
  QCheckBox* WindowLevelCheckBox;
  
  /// The smallest value in the image.
  float RangeFrom;

  /// The largest value in the image.
  float RangeTo;

  /// The dominant width of the range (standard deviation or total width).
  float RangeWidth;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtWindowLevelControls_h_included_

