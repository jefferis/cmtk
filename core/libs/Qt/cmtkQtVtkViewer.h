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

#include <QGridLayout>

#ifndef __cmtkQtVtkViewer_h_included_
#define __cmtkQtVtkViewer_h_included_

#include <cmtkconfig.h>

#ifdef CMTK_HAVE_VTK

#include <qwidget.h>
#include <qslider.h>
#include <qlayout.h>
#include <qmenubar.h>

#include <vtkQtRenderWindow.h>
#include <vtkQtRenderWindowInteractor.h>

#include <vtkRenderer.h>
#include <vtkWindowToImageFilter.h>
#include <vtkTIFFWriter.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

/// Three-dimensional viewer for fusion application.
class QtVtkViewer :
  public QWidget
{
  Q_OBJECT

public:
  /// Constructor.
  QtVtkViewer
  ( QWidget *const parent = 0, const char* name = 0, Qt::WFlags flags = 0 );

  /// Constructor with explicit numbers of layout rows and columns.
  QtVtkViewer
  ( const int nrows, const int ncolumns, 
    QWidget *const parent = 0, const char* name = 0, Qt::WFlags flags = 0 );

  /// Destructor.
  virtual ~QtVtkViewer() {}

  /// Return pointer to renderer.
  vtkRenderer* GetRenderer( const int row = 1, const int col = 1 ) {
    return Renderer[row][col];
  }

public slots:
  /// Reset renderer camera.
  void slotResetCamera();

  /// Reset renderer camera to AP direction.
  void slotSetCameraDirectionAP();

  /// Reset renderer camera to AP direction.
  void slotSetCameraDirectionLR();

  /// Reset renderer camera to AP direction.
  void slotSetCameraDirectionHF();

  /// Reset renderer camera to AP direction.
  void slotSetCameraDirectionPA();

  /// Reset renderer camera to AP direction.
  void slotSetCameraDirectionRL();

  /// Reset renderer camera to AP direction.
  void slotSetCameraDirectionFH();

  /// Set camera azimuth.
  void slotAzimuth( int angle );

  /// Set camera elevation.
  void slotElevation( int angle );

  /// Set camera dolly.
  void slotDolly( float dolly );

  /// Increase dolly by 10.
  void slotDollyPlus10() { this->slotDolly( 1.1 ); }
  /// Decrease dolly by 10.
  void slotDollyMinus10() { this->slotDolly( 0.9 ); }

  /// Increase azimuth by 10.
  void slotAzimuthPlus10() { this->slotAzimuth( 10 ); }
  /// Decrease azimuth by 10.
  void slotAzimuthMinus10() { this->slotAzimuth( -10 ); }

  /// Increase elevation by 10.
  void slotElevationPlus10() { this->slotElevation( 10 ); }
  /// Decrease elevation by 10.
  void slotElevationMinus10() { this->slotElevation( -10 ); }

  /// Export render view to image file.
  void slotExportImage( const QString& filename );

protected slots:
  /// Update all visible renderers.
void slotUpdateRenderers();

private slots:
  /// Set window layout to single viewer.
  void slotSetLayoutSingle();

  /// Set window layout to panel of four viewers.
  void slotSetLayoutPanel();

public:
  /// Add actor to scene.
  void AddActor( vtkActor *const actor, const bool render = true ) 
  {
    for ( int row = 0; row < 2; ++row )
      for ( int col = 0; col < 2; ++col )
	Renderer[row][col]->AddActor( actor );
    if ( render ) 
      this->slotUpdateRenderers();
  }

protected:
  /// The top-level layout.
  QGridLayout* MasterLayout;

  /// Layout of render windows.
  QGridLayout* RenderLayout;

  /// Menu bar.
  QMenuBar* MenuBar;

  /// The VTK renderer.
  vtkRenderer* Renderer[2][2];

  /// Render window.
  vtkQtRenderWindow* RenderWindow;

  /// Render window interactor.
  vtkQtRenderWindowInteractor* Interactor;

  /// Set camera based on view normal and view up vector.
  void SetCamera( const double* vn, const double* vup, const int row = 1, const int col = 1 );

private:
  // Construct dialog with given numbers of rows and columns in layout.
  void Construct( const int nrows = 2, const int ncols = 3 );

  /// Window-to-image filter.
  vtkWindowToImageFilter* WindowToImage;

  /// VTK TIFF writer.
  vtkTIFFWriter* TIFFWriter;
};

//@}

} // namespace cmtk

#endif // #ifdef CMTK_HAVE_VTK

#endif // #ifndef __cmtkQtVtkViewer_h_included_
