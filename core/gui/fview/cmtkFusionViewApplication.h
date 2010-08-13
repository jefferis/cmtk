/*
//
//  Copyright 2010 SRI International
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

#ifndef __cmtkFusionViewApplication_h_included_
#define __cmtkFusionViewApplication_h_included_

#include <QtGui/QApplication>
#include <QtGui/QMainWindow>

#include <Base/cmtkUniformVolume.h>

#include <ui_fviewMainWindow.h>

namespace
cmtk
{

/// Application class for fusion viewer.
class FusionViewApplication : public QApplication
{
  Q_OBJECT

public:
  /// Constructor.
  FusionViewApplication( int argc, char* argv[] );

private:
  /// Application main window.
  QMainWindow* m_MainWindow;

  /// Designed-generated User Interface for the main window.
  Ui::fviewMainWindow m_MainWindowUI;

  /// The fixed volume.
  UniformVolume::SmartConstPtr m_FixedVolume;

  /// The moving volume.
  UniformVolume::SmartConstPtr m_MovingVolume;
};

} // namespace cmtk

#endif // #ifndef __cmtkFusionViewApplication_h_included_
