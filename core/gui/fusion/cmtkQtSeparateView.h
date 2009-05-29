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

#ifndef QtSeparateView_h_included_
#define QtSeparateView_h_included_

#include <cmtkconfig.h>

#include <cmtkQtFusionWindowTemplate.h>

#include <cmtkQtSimpleFusionApp.h>
#include <cmtkStudyList.h>

#include <cmtkQtScrollRenderView.h>
#include <cmtkImageToImageRGB.h>
#include <cmtkColormap.h>

#include <cmtkQtStudyNamesBox.h>

namespace
cmtk
{

/// Class for separate display of registered images.
class QtSeparateView :
  /// Inherit from fusion window template.
  public QtFusionWindowTemplate
{
  Q_OBJECT

public:
  /// Constructor.
  QtSeparateView( QtSimpleFusionApp* fusionApp, QWidget *const parent = 0, Qt::WFlags flags = 0 );
  
  /** Destructor.
   */
  virtual ~QtSeparateView() {};

  /// Export image: This function needs to be overwritten by derived classes.
  virtual void Export( const QString&, const QString& = QString::null, const QStringList* = NULL ) {};

public slots:
  /// React to changed slice plane.
  void slotUpdateSlice();

  /// React to changed study colormap.
  void slotUpdateColormap( Study::SmartPtr& study );

private slots:
  /// Update study colormaps.
  void slotUpdateColormaps();

  /// Switch displayed study in left window.
  void slotSwitchStudyL( const QString& studyName );

  /// Switch displayed study in right window.
  void slotSwitchStudyR( const QString& studyName );

protected:
  /// Update list widgets with study names.
  virtual void UpdateStudySelection();

  /// This virtual member is called when the slice changes.
  virtual void UpdateSlice() 
  { 
    Left.View->slotRender(); 
    Right.View->slotRender(); 
  }

private:
  class UI 
  {
  public:
    /// Constructor.
    UI() : m_Study( NULL ) {}

    /// Destructor.
    ~UI() 
    { 
      this->m_ImageToImageRGB->SetInput( NULL ); 
      this->m_ImageToImageRGB->Delete();
      this->m_Colormap->Delete(); 
    }
    
    /// Construct this object.
    void Construct( QWidget *const parent, QLayout *const inLayout, const QString& label );

    /// Link to the study object.
    Study::SmartPtr m_Study;

    /// The view widget.
    QtScrollRenderView* View;

    /// The current colormap.
    Colormap* m_Colormap;
    
    /// The colormap lookup filter.
    ImageToImageRGB* m_ImageToImageRGB;    

    /// The box widget that contains all available study names.
    QtStudyNamesBox* StudyNamesBox;
  };
  
  UI Left;
  UI Right;
};

} // namespace cmtk

#endif // #ifndef QtSeparateView_h_included_
