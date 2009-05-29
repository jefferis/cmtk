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

#ifndef QtFusionAlpha_h_included_
#define QtFusionAlpha_h_included_

#include <cmtkconfig.h>

#include <cmtkQtFusionWindowTemplate.h>
#include <cmtkQtSimpleFusionApp.h>
#include <cmtkQtScrollRenderView.h>

#include <qslider.h>

#include <cmtkStudyList.h>

#include <cmtkImageToImageRGB.h>
#include <cmtkColormap.h>
#include <cmtkFusionAlpha.h>

#include <cmtkQtStudyNamesBox.h>

namespace
cmtk
{

/// Class for alpha blending image fusion.
class QtFusionAlpha :
  /// Inherit from fusion window template.
  public QtFusionWindowTemplate
{
  Q_OBJECT

public:
  /// Constructor.
  QtFusionAlpha( QtSimpleFusionApp *const fusionApp, QWidget *const parent = 0, Qt::WFlags flags = 0 );

  /** Destructor.
   * Takes care of destroying the internal visualization pipelines.
   */
  virtual ~QtFusionAlpha();

  /// Flag whether top study is also alpha (transparency) study.
  bool TopIsAlpha;

public slots:
  /// React to changed slice plane.
  void slotUpdateSlice();
  
  /// React to changed study colormap.
  void slotUpdateColormap( Study::SmartPtr& study );

  /// React to "Update" button.
  void slotUpdateColormaps();

private slots:
  /// Switch top study.
  void slotSwitchStudyTop( const QString& studyName );

  /// Switch bottom study.
  void slotSwitchStudyBottom( const QString& studyName );

  /// Switch alpha study.
  void slotSwitchStudyAlpha( const QString& studyName );

  /// Value of "ramp from" slider changed.
  void slotSliderFromChanged( double value );

  /// Value of "ramp to" slider changed.
  void slotSliderToChanged( double value );

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

  Study::SmartPtr StudyTop;
  Study::SmartPtr StudyBottom;
  Study::SmartPtr StudyAlpha;

  Colormap* ColormapTop;
  Colormap* ColormapBottom;
  Colormap* ColormapAlpha;

  ImageToImageRGB* ImageToImageRGBTop;
  ImageToImageRGB* ImageToImageRGBBottom;
  ImageToImageRGB* ImageToImageRGBAlpha;

  FusionAlpha* FusionFilter;

  /// Top list box with available study names.
  QtStudyNamesBox* StudyNamesBoxTop;

  /// Bottom list box with available study names.
  QtStudyNamesBox* StudyNamesBoxBottom;

  /// Alpha list box with available study names.
  QtStudyNamesBox* StudyNamesBoxAlpha;

  QtSliderEntry* SliderFrom;
  QtSliderEntry* SliderTo;
};

} // namespace cmtk

#endif // #ifndef QtFusionAlpha_h_included_
