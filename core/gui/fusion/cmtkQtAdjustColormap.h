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

#ifndef __cmtkQtAdjustColormap_h_include_
#define __cmtkQtAdjustColormap_h_include_

#include <cmtkconfig.h>

#include <qwidget.h>

#include <cmtkMacros.h>
#include <cmtkStudy.h>

#include <cmtkImage.h>
#include <cmtkColormap.h>
#include <cmtkImageToImageRGB.h>

#include <cmtkQtScrollRenderView.h>
#include <cmtkQtWindowLevelControls.h>

#include <QEvent>
#include <QPixmap>

namespace
cmtk
{

/** Image browser window.
 */
class QtAdjustColormap :
  public QWidget
{
  Q_OBJECT // we're using slots

public:
  /// Constructor.
  QtAdjustColormap();

  /// Destructor.
  virtual ~QtAdjustColormap();

  /// Switch to another study.
  virtual void SwitchToStudy( Study::SmartPtr& study );

private:
  /// The Study object we're working on.
  cmtkGetSetMacro(Study::SmartPtr,Study);

  /// Update dialog after study change.
  void UpdateDialog();

signals:
  /// This signal is emitted whenever the colormap changes.
  void colormapChanged( const cmtk::Study* study );

  /// This signal is emitted when a new icon was generated for a study.
  void newIcon( const Study* study, QPixmap icon );

private slots:
  /// Switch image.
  void slotSwitchImage( int imageIndex );

  /// This slot handles changes made to the window/level UI.
  void slotWindowChanged();

  /// This slot notifies the browser of a palette selection event.
  void slotSetStandardColormap( int paletteIdx );

protected:
  /// Handle "focus lost" events.
  virtual void leaveEvent( QEvent *event );

private:
  /// The scrolled view we display an image in.
  QtScrollRenderView* ScrollRenderView;

  QtWindowLevelControls* WindowLevelBox;

  Image* PipelineImage;

  Colormap* m_Colormap;

  ImageToImageRGB* m_ImageToImageRGB;
};

} // namespace cmtk

#endif // #ifndef __cmtkQtAdjustColormap_h_include_
