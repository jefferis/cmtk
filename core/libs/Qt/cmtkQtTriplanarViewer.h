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

#ifndef __cmtkQtTriplanarViewer_h_included_
#define __cmtkQtTriplanarViewer_h_included_

#include <cmtkconfig.h>

#include <cmtkQtTriplanarWindow.h>
#include <cmtkQtWindowLevelDialog.h>

#include <qlistwidget.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/** Stand-alone triplanar image viewer.
 */
class QtTriplanarViewer :
  /// Inherit from triplanar viewer widget.
  public QtTriplanarWindow
{
  Q_OBJECT // we're using slots

public:
  /// Constructor.
  QtTriplanarViewer();

  /// Virtual destructor.
  virtual ~QtTriplanarViewer() {};

  /// Execute in batch mode.
  virtual int ExecuteBatchMode( const int argc, char* argv[] );

public slots:
  /// Add study by filesystem path.
  void slotAddStudy( const char* fname );

  /// Load image from file.
  void slotLoadFile();

  /// Load image from file.
  void slotReloadData();

  /// Copy current image colormap to all other images.
  void slotCopyColormapToOtherImages();

private:
  /// Window/Level dialog.
  QtWindowLevelDialog* WindowLevel;

  /// Vector of loaded studies.
  std::vector<Study::SmartPtr> m_Studies;

  /// Tab for the images list.
  QWidget* m_ImagesTab;

  /// List box with loaded studies' names.
  QListWidget* m_StudiesListBox;

private slots:
  /// Study was double-clicked in listbox.
  void slotSwitchStudy( const QString & study );
};

//@}

} // namespace cmtk

#endif // #ifndef __cmtkQtTriplanarViewer_h_included_
