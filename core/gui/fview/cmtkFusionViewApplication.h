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

  /// The slice axis (0=x, sagittal; 1=y, coronal; 2=z, axial).
  int m_SliceAxis;

  /// Slice index in the fixed image along the slice axis.
  int m_SliceIndex;

  /// Data for the current fixed image slice.
  UniformVolume::SmartConstPtr m_FixedSlice;

  /// Update displayed fixed image slice.
  void UpdateFixedSlice();

  /// Update displayed moving image slice.
  void UpdateMovingSlice();

  /// Update widget using slice data, black and white levels.
  void UpdateWidget( QWidget* widget, const UniformVolume& slice, const float blackLevel, const float whiteLevel );
};

} // namespace cmtk

#endif // #ifndef __cmtkFusionViewApplication_h_included_
