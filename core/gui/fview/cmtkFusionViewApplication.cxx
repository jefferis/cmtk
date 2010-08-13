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

#include "cmtkFusionViewApplication.h"

#include <ui_fviewMainWindow.h>

#include <System/cmtkCommandLine.h>

#include <IO/cmtkVolumeIO.h>

cmtk::FusionViewApplication
::FusionViewApplication( int argc, char* argv[] ) 
  : QApplication( argc, argv ),
    m_MainWindow( new QMainWindow )
{
  CommandLine cl;
  cl.SetProgramInfo( CommandLine::PRG_TITLE, "Fusion viewer." );

  const char* imagePathFix;
  cl.AddParameter( &imagePathFix, "FixedImage", "Fixed image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

  const char* imagePathMov;
  cl.AddParameter( &imagePathMov, "MovingImage", "Moving image path" )->SetProperties( cmtk::CommandLine::PROPS_IMAGE );

  try
    {
    cl.Parse( argc, const_cast<const char**>( argv ) );
    }
  catch ( const CommandLine::Exception& ex )
    {
    throw(ex);
    }

  this->m_FixedVolume = VolumeIO::Read( imagePathFix );
  if ( ! this->m_FixedVolume )
    {
    exit( 1 );
    }

  this->m_MovingVolume = VolumeIO::Read( imagePathMov );
  if ( ! this->m_MovingVolume )
    {
    exit( 1 );
    }
  
  Ui::fviewMainWindow ui;
  ui.setupUi( this->m_MainWindow );
  this->m_MainWindow->show();
}

