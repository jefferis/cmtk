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

#ifndef __cmtkQtSimpleFusionApp_h_included_
#define __cmtkQtSimpleFusionApp_h_included_

#include <cmtkconfig.h>

#include <qapplication.h>

#include <cmtkQtSimpleFusionMainWindow.h>

#include <cmtkResourceFile.h>
#include <cmtkStudy.h>
#include <cmtkStudyList.h>
#include <cmtkQtFusionSlicer.h>
#include <QLabel>

namespace
cmtk
{

/** Image fusion application class.
 */
class QtSimpleFusionApp :
  /// Inherit from Qt's application class.
  public QApplication
{
  Q_OBJECT

public:
  /// The current reference study.
  Study::SmartPtr ReferenceStudy;

  /// The current study list.
  StudyList::SmartPtr m_StudyList;

public slots:
  /// Slot for setting new studylist.
  void slotSetStudyList( StudyList::SmartPtr& studyList );

  /// Slot for adding a new study.
  void slotAddStudy( Study::SmartPtr& study );

  /// Slot for adding a new study.
  void slotSetReferenceStudy( Study::SmartPtr& study );

  /// Slot for notification when image data changes.
  void slotDataChanged( Study::SmartPtr& study );

signals:
  /// Slice changed.
  void signalSliceChanged();

  /// Studylist changed.
  void signalStudyListChanged();

  /// Reference study changed.
  void signalReferenceStudyChanged();

  /// Image data was changed.
  void signalDataChanged( Study::SmartPtr& study );

public:
  /// Constructor.
  QtSimpleFusionApp( int argc, char *argv[] );

  /// Destructor.
  virtual ~QtSimpleFusionApp();

  /// Resource file. 
  ResourceFile m_ResourceFile;

  /// Get (potentially create) fusion slicer.
  QtFusionSlicer* GetFusionSlicer();

private:
  /// The actual resource file name including path.
  QString ResourceFilePath;
  
  /// The application's main window.
  QtSimpleFusionMainWindow* MainWindow;

  /// Fusion slicer object.
  QtFusionSlicer* m_FusionSlicer;
};

} // namespace cmtk

#endif // #ifndef __cmtkQtSimpleFusionApp_h_included_
