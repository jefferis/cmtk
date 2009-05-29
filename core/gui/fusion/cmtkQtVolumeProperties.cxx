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

#include <cmtkconfig.h>

#include <cmtkQtVolumeProperties.h>

#include <qlayout.h>
#include <qstring.h>
#include <qlabel.h>
#include <q3groupbox.h>
#include <qpushbutton.h>
#include <Q3VBoxLayout>

#include <cmtkQtFusionGlobal.h>

namespace
cmtk
{

QtVolumeProperties::QtVolumeProperties( const Study* study )
  : QWidget( NULL, "VolumeProperties", 0 )
{
  this->setIcon( QtFusionGlobal::WindowIcon() );
  Q3VBoxLayout* layout = new Q3VBoxLayout( this );
    
  if ( study ) 
    {    
    QString caption;
    this->setCaption( caption.sprintf( "Image Properties: %s", study->GetName() ) );

    Q3GroupBox* generalBox = new Q3GroupBox( 2, Qt::Horizontal, this );
    generalBox->setTitle( "General" );
    layout->addWidget( generalBox );
    
    new QLabel( "Path: ", generalBox );
    new QLabel( study->GetFileSystemPath(), generalBox );

    new QLabel( "Description: ", generalBox );
    new QLabel( study->GetDescription(), generalBox );

    Q3GroupBox* gridBox = new Q3GroupBox( 2, Qt::Horizontal, this );
    gridBox->setTitle( "Coordinate Grid" );
    layout->addWidget( gridBox );

    const Volume *volume = study->GetVolume();
    if ( volume ) 
      {
      const UniformVolume *uniVolume = dynamic_cast<const UniformVolume*>( volume );
      if ( uniVolume ) 
	{
	new QLabel( "Type: ", gridBox );
	new QLabel( "Uniform", gridBox );

	new QLabel( "Dimensions: ", gridBox );
	QString dims;
	dims.sprintf( "%d x %d x %d pixels", uniVolume->GetDims(AXIS_X), uniVolume->GetDims(AXIS_Y), uniVolume->GetDims(AXIS_Z) );
	new QLabel( dims, gridBox );

	new QLabel( "Total size: ", gridBox );
	QString size;
	size.sprintf( "%f x %f x %f mm", uniVolume->Size[AXIS_X], uniVolume->Size[AXIS_Y], uniVolume->Size[AXIS_Z] );
	new QLabel( size, gridBox );

	new QLabel( "Pixel size: ", gridBox );
	QString psize;
	psize.sprintf( "%f x %f x %f mm", uniVolume->Delta[AXIS_X], uniVolume->Delta[AXIS_Y], uniVolume->Delta[AXIS_Z] );

	new QLabel( psize, gridBox );
	} 
      else
	{
	new QLabel( "Type: ", gridBox );
	new QLabel( "Non-uniform rectilinear", gridBox );	
	}
      }
    }
  
  QPushButton* closeButton = new QPushButton( "Close", this );
  layout->addWidget( closeButton );
  QObject::connect( closeButton, SIGNAL( clicked() ), this, SLOT( close() ) );

  this->show();
}

void QtVolumeProperties::slotReadData()
{
}

} // namespace cmtk
