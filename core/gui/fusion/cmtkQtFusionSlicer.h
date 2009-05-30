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

#ifndef __cmtkQtFusionSlicer_h_included_
#define __cmtkQtFusionSlicer_h_included_

#include <cmtkconfig.h>

#include <qwidget.h>
#include <qcombobox.h>
#include <qstringlist.h>

#include <cmtkTypes.h>
#include <cmtkInterpolator.h>

#include <cmtkStudy.h>
#include <cmtkStudyList.h>
#include <cmtkFusionSlicers.h>

#include <cmtkQtStudyNamesBox.h>
#include <cmtkQtSliderEntry.h>

namespace
cmtk
{

class QtSimpleFusionApp;

/// User interface for planar slicer.
class QtFusionSlicer :
  /// This is a generic window widget.
  public QWidget
{
  Q_OBJECT

public:
  /// Constructor.
  QtFusionSlicer( QtSimpleFusionApp *const fusionApp );

  /// Virtual destructor.
  ~QtFusionSlicer();

public slots:
  /// Set studylist.
  void slotStudyListChanged();

  /// Set reference study.
  void slotReferenceStudyChanged();

  /// Set reference study.
  void slotSetReferenceStudy( const QString& studyName );

  /// Set plane orientation.
  void slotSetOrientation( int sliceNormal );

  /// Set slice position.
  void slotSetSlicePosition( double slicePosition );

  /// Set interpolation mode.
  void slotSetInterpolation( int mode );

  /// Set interpolation mode.
  void slotApplyWarp( int mode );

signals:
  /// This signal is emitted when the widget sets a new reference study.
  void signalNewReferenceStudy( Study::SmartPtr& );

  /// This signal is emitted when the slice plane has changed.
  void sliceChanged();

  /** This signals all viewer windows to write their contents to a file.
   *\param path File system path for the output file.
   *\param format Image file format.
   @\param comments List of file comments.
  */
  void signalExport( const QString&, const QString&, const QStringList* );

public:
  /// Get a constant pointer to the list of study names.
  const QStringList* GetStudyNamesList() const 
  {
    return StudyNamesList;
  }

  /// Get a constant pointer to the list of target study names.
  const QStringList* GetTargetStudyNamesList() const 
  {
    return TargetStudyNamesList;
  }

  /// Get given output slice plane.
  Image* GetOutput( Study::SmartPtr& study );

private:
  /// Link to the fusion application.
  QtSimpleFusionApp* FusionApp;

  /// This is the low-level object that handles the actual slice pipelines.
  FusionSlicers* m_FusionSlicers;

  /// List of all available study names for convenient access.
  QStringList* StudyNamesList;

  /// List of names of studies registered to current reference study.
  QStringList* TargetStudyNamesList;

  /// Combobox for reference study selection.
  QtStudyNamesBox* ReferenceBox;

  /// Slider for plane position.
  QtSliderEntry* SliceSlider;

  /// Slice normal (axis perpendicular to slice plane).
  int SliceNormal;

  /// Maximum slice position.
  Types::Coordinate PositionTo;

  /// Minimal increment between two slices.
  Types::Coordinate PositionStep;
};

} // namespace cmtk

#endif // #ifndef __cmtkQtFusionSlicer_h_included_
