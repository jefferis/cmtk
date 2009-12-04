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

#ifndef __cmtkQtImageOperators_h_included_
#define __cmtkQtImageOperators_h_included_

#include <cmtkconfig.h>

#include <qobject.h>
#include <qwidget.h>
#include <qmenu.h>

#include <cmtkStudy.h>
#include <cmtkQtProgress.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/// A collection of 3D image operators with Qt menu.
class QtImageOperators :
  /// Inherit from QObject for event handling etc.
  public QObject
{
  Q_OBJECT

public slots:
  /// Slot to open selected type of fusion window.
  void slotOperatorsMedian();
  void slotOperatorsSobel();
  void slotOperatorsHistEq();
  void slotOperatorsAbs();
  void slotOperatorsLog();
  void slotOperatorsExp();

signals:
  /// This signal is sent when the image data has been changed.
 void dataChanged( Study::SmartPtr& );

public:
  /// Constructor.
  QtImageOperators
  ( Study::SmartPtr* currentStudy, QWidget *const parent = NULL, QtProgress *const progressInstance = NULL )
    : Parent( parent ), CurrentStudy( currentStudy ),
      ProgressInstance( progressInstance ) {};
  
  /// Create and return popup menu that makes operators available.
  QMenu* CreatePopupMenu();

private:
  /** The parent widget.
   * This is for modal dialogs that may be opened for some operations.
   */
  QWidget* Parent;

  /// Pointer to an object with the current study pointer.
  Study::SmartPtr* CurrentStudy;

  /// Optional instance of a Qt progress indicator.
  QtProgress* ProgressInstance;

  /// Check whether study and volume data are all valid.
  bool StudyDataValid() const 
  {
    return (*(this->CurrentStudy)) && (*(this->CurrentStudy))->GetVolume() && (*(this->CurrentStudy))->GetVolume()->GetData();
  }
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtImageOperators_h_included_
