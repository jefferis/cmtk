/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#ifndef QtFusionEdge_h_included_
#define QtFusionEdge_h_included_

#include <cmtkconfig.h>

#include "cmtkQtFusionWindowTemplate.h"
#include "cmtkQtSimpleFusionApp.h"
#include "cmtkQtStudyNamesBox.h"

#include "Qt/cmtkQtScrollRenderView.h"
#include "Qt/cmtkQtWindowLevelControls.h"
#include "Qt/cmtkQtSliderEntry.h"

#include "IO/cmtkStudyList.h"

#include "Pipeline/cmtkImageToImageRGB.h"
#include "Pipeline/cmtkColormap.h"
#include "Pipeline/cmtkFusionEdge.h"

#include <qslider.h>
#include <qgroupbox.h>
#include <qcheckbox.h>


namespace
cmtk
{

/// Class for alpha blending image fusion.
class QtFusionEdge :
  /// Inherit from fusion window template.
  public QtFusionWindowTemplate
{
  Q_OBJECT

public:
  /// Constructor.
  QtFusionEdge( QtSimpleFusionApp *const fusionApp, QWidget *const parent = 0, Qt::WFlags flags = 0 );

  /** Destructor.
   * Takes care of destroying the internal visualization pipelines.
   */
  virtual ~QtFusionEdge();

public slots:
  /// React to changed slice plane.
  void slotUpdateSlice();
  
  /// React to changed study colormap from the outside.
  void slotUpdateColormap( Study::SmartPtr& study );

  /// Update colormaps from global objects ("Update" button pressed).
  void slotUpdateColormaps();

private slots:
  /// Switch top study.
  void slotSwitchStudyEdge( const QString& studyName );

  /// Switch bottom study.
  void slotSwitchStudyBG( const QString& studyName );

  /// Value of "ramp from" slider changed.
  void slotSliderFromChanged( int value );

  /// Value of "ramp to" slider changed.
  void slotSliderToChanged( int value );

  /// Different edge colormap selected.
  void slotSetColormap( int value );

  /// Parameters in the UI changed.
  void slotParametersChanged();

protected:
  /// Update list widgets with study names.
  virtual void UpdateStudySelection();

  /// This virtual member is called when the slice changes.
  virtual void UpdateSlice() { View->slotRender(); };

  /// Export displayed image.
  virtual void Export
  ( const QString& path, const QString& format = QString::null, 
    const QStringList* comments = NULL );

private:
  QtScrollRenderView* View;
  QtWindowLevelControls* WindowLevelBox;
  Study::SmartPtr StudyEdge;
  Study::SmartPtr StudyBG;
  Study::SmartPtr StudyEdgeImage;

  FusionEdge* m_FusionEdge;

  /// Top list box with available study names.
  QtStudyNamesBox* StudyNamesBoxEdge;

  /// Bottom list box with available study names.
  QtStudyNamesBox* StudyNamesBoxBG;

  QSlider* SliderFrom;
  QSlider* SliderTo;
  QGroupBox* OperatorButtons;
  QComboBox* OperatorBox;
  QCheckBox* SmoothCheckBox;
  QtSliderEntry* GaussianWidthSlider;

  int EdgeStandardColormap;
};

} // namespace cmtk

#endif // #ifndef QtFusionEdge_h_included_
