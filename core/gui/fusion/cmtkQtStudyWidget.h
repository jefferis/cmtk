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

#ifndef __cmtkQtStudyWidget_h_included_
#define __cmtkQtStudyWidget_h_included_

#include <cmtkconfig.h>

#include <qwidget.h>

#include <cmtkStudy.h>

#include <cmtkImage.h>
#include <cmtkColormap.h>
#include <cmtkImageToImageRGB.h>

#include <cmtkQtScrollRenderView.h>
#include <cmtkQtWindowLevelControls.h>

namespace
cmtk
{

/// Study controls widget, for example for the main window tabs.
class QtStudyWidget :
  /// Inherit generic Qt widget.
public QWidget
{
  Q_OBJECT

public:
  /// Constructor
  QtStudyWidget( QWidget* parent, const char* name = 0, Qt::WFlags flags = 0 );

  /// Get study object.
  Study::SmartPtr& GetStudy() 
  { 
    return this->m_Study; 
  }

signals:
  void colormap( Study::SmartPtr& study );
  void volume( Study::SmartPtr& study );
  void deleted( Study::SmartPtr& study );

public slots:
  /// Slot to set study.
  void slotSetStudy( Study::SmartPtr& study );

  /// Slot update view after image data changed.
  void slotDataChanged( Study::SmartPtr& study );

private slots:
  /// Slot to switch to other image.
  void slotSwitchImage( int imageIndex );

  /// Slot called when window/level sliders change.
  void slotUpdateColormap();

private:
  /// Smart pointer to study object.
  Study::SmartPtr m_Study;

  /// The scrolled view we display an image in.
  QtScrollRenderView* ScrollRenderView;

  QtWindowLevelControls* WindowLevelBox;

  Image* PipelineImage;

  Colormap* m_Colormap;

  ImageToImageRGB* m_ImageToImageRGB;

  /// Index of the currently displayed slice.
  int ImageIndex;
};

} // namespace cmtk

#endif // #ifndef __cmtkQtStudyWidget_h_included_
