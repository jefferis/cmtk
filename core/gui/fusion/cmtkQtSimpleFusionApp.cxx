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

#include <cmtkQtSimpleFusionApp.h>

#include <string.h>

#include <qpixmap.h>
#include <qimage.h>
#include <qstring.h>
#include <qdir.h>

#include <QLabel>
#include <Q3Frame>

namespace
cmtk
{

QtSimpleFusionApp::QtSimpleFusionApp( int argc, char *argv[] )
  : QApplication( argc, argv ),
    ReferenceStudy( NULL ), 
    m_StudyList( new StudyList ),
    m_FusionSlicer( NULL )
{
  QString rcFileName = ".fusionrc";
  QDir qdir = QDir::current();
  if ( qdir.exists( rcFileName ) ) 
    {
    rcFileName = qdir.absFilePath( rcFileName );
    } 
  else
    {
    qdir = QDir::home();
    rcFileName = qdir.absFilePath( rcFileName );
    }
  ResourceFilePath = rcFileName;
  this->m_ResourceFile.Read( ResourceFilePath.latin1() );
  
  MainWindow = new QtSimpleFusionMainWindow( this );
  this->setMainWidget( MainWindow );
  
  ResourceSection& section = this->m_ResourceFile["MainWindow"];
  ResourceSection::const_iterator it = section.begin();

  while ( it != section.end() ) 
    {
    if ( !strcmp( it->c_str(), "maximized = yes" ) ) 
      {
      MainWindow->showMaximized();
      }
    
    int x, y, w, h;
    if ( 4 == sscanf( it->c_str(), "geometry %d %d %d %d", &x, &y, &w, &h ) ) 
      {
      MainWindow->setGeometry( x, y, w, h );
      }
    ++it;
    }
  
  MainWindow->show();
  ///  delete splash;
}

QtSimpleFusionApp::~QtSimpleFusionApp() 
{
  ResourceSection& section = this->m_ResourceFile["MainWindow"];
  section.clear();

  if ( MainWindow->isMaximized() )
    section.push_back( "maximized = yes" );
  else
    section.push_back( "maximized = no" );
          
  QString geomStr;
  geomStr.sprintf( "geometry %d %d %d %d", MainWindow->x(), MainWindow->y(), MainWindow->width(), MainWindow->height() );
  section.push_back( geomStr.latin1() );

  this->m_ResourceFile.Write( ResourceFilePath );
}

QtFusionSlicer* 
QtSimpleFusionApp::GetFusionSlicer()
{
  if ( ! this->m_FusionSlicer ) 
    {
    this->m_FusionSlicer = new QtFusionSlicer( this );
    }
  return this->m_FusionSlicer;
}

} // namespace cmtk
