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

#include <cmtkQtTriplanarViewer.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

int
QtTriplanarViewer
::ExecuteBatchMode( const int argc, char* argv[] )
{
  this->m_BatchMode = true;
  for ( int i = 0; i < argc; ++i )
    {
    if ( !strcmp( argv[i], "load" ) )
      {
      this->slotAddStudy( argv[++i] );      
      }
    else if ( !strcmp( argv[i], "goto-pixel" ) )
      {
      this->slotGoToPixel( argv[++i] );
      }
    else if ( !strcmp( argv[i], "goto-location" ) )
      {
      this->slotGoToLocation( argv[++i] );
      }
    else if ( !strcmp( argv[i], "colormap" ) )
      {
      this->slotSetColormap( argv[++i] );
      }
    else if ( !strcmp( argv[i], "window-level" ) )
      {
      this->slotSetWindowLevel( argv[++i] );
      }
    else if ( !strcmp( argv[i], "zoom" ) )
      {
      this->slotSetZoom( atoi( argv[++i] ) );
      }
    else if ( !strcmp( argv[i], "crosshair" ) )
      {
      const char* chOnOff = argv[++i];
      this->slotSetCrosshairMode( ! strcmp( chOnOff, "on" ) || ! strcmp( chOnOff, "yes" ) || ! strcmp( chOnOff, "true" ) );
      }
    else if ( !strcmp( argv[i], "checkerboard" ) )
      {
      const char* chOnOff = argv[++i];
      this->slotSetCheckerboardMode( ! strcmp( chOnOff, "on" ) || ! strcmp( chOnOff, "yes" ) || ! strcmp( chOnOff, "true" ) );
      }
    else if ( !strcmp( argv[i], "export-axial" ) )
      {
      this->slotExportImage( argv[++i], 1 );
      }
    else if ( !strcmp( argv[i], "export-coronal" ) )
      {
      this->slotExportImage( argv[++i], 2 );
      }
    else if ( !strcmp( argv[i], "export-sagittal" ) )
      {
      this->slotExportImage( argv[++i], 3 );
      }
    else if ( !strcmp( argv[i], "export-panel" ) )
      {
      this->slotExportImage( argv[++i], 4 );
      }
    }
  return 0;
}

} // namespace cmtk
