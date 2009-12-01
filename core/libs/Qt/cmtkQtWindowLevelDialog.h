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

#ifndef __cmtkQtWindowLevelDialog_h_included_
#define __cmtkQtWindowLevelDialog_h_included_

#include <cmtkconfig.h>

#include <qdialog.h>

#include <cmtkQtWindowLevelControls.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/// Dialog with WIndow/Level and Colormap controls.
class QtWindowLevelDialog :
  /// Inherit from Qt dialog.
  public QDialog
{
  Q_OBJECT

public:
  /// Constructor.
  QtWindowLevelDialog( QWidget* parent = 0, const char* name = 0, bool modal = FALSE, Qt::WFlags f = 0 );

  /// Virtual destructor.
  virtual ~QtWindowLevelDialog() {}

public slots:
  /// Set study object.
  void slotSetStudy( Study::SmartPtr& study );

signals:
  /// This signal is emitted when the colormap of the study has changed.
  void colormapChanged( Study::SmartPtr& );

private:
  /// The Window/Level user interface component.
  QtWindowLevelControls* Controls;
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtWindowLevelDialog_h_included_
