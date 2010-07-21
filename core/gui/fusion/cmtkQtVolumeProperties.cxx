/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include <cmtkconfig.h>

#include "cmtkQtVolumeProperties.h"
#include "cmtkQtFusionGlobal.h"

#include <qlayout.h>
#include <qstring.h>
#include <qlabel.h>
#include <qgroupbox.h>
#include <qpushbutton.h>
#include <QGridLayout>

namespace
cmtk
{

QtVolumeProperties::QtVolumeProperties( const Study* study )
  : QWidget( NULL )
{
  this->setWindowIcon( QtFusionGlobal::WindowIcon() );
  QVBoxLayout* layout = new QVBoxLayout( this );
    
  if ( study ) 
    {    
    QString caption;
    this->setWindowTitle( caption.sprintf( "Image Properties: %s", study->GetName() ) );

    QGroupBox* generalBox = new QGroupBox( this );
    QGridLayout *generalLayout = new QGridLayout;

    generalBox->setTitle( "General" );
    generalBox->setLayout( generalLayout );
    layout->addWidget( generalBox );
    
    generalLayout->addWidget( new QLabel( "Path: ", generalBox ), 0, 0 );
    generalLayout->addWidget( new QLabel( study->GetFileSystemPath(), generalBox ), 0, 1 );

    generalLayout->addWidget( new QLabel( "Description: ", generalBox ), 1, 0 );
    generalLayout->addWidget( new QLabel( study->GetDescription(), generalBox ), 1, 1 );

    QGroupBox* gridBox = new QGroupBox( this );
    QGridLayout *gridLayout = new QGridLayout;

    gridBox->setTitle( "Coordinate Grid" );
    gridBox->setLayout( gridLayout );
    layout->addWidget( gridBox );

    const Volume *volume = study->GetVolume();
    if ( volume ) 
      {
      const UniformVolume *uniVolume = dynamic_cast<const UniformVolume*>( volume );
      if ( uniVolume ) 
	{
	gridLayout->addWidget( new QLabel( "Type: ", gridBox ), 0, 0 );
	gridLayout->addWidget( new QLabel( "Uniform", gridBox ), 0, 1 );

	gridLayout->addWidget( new QLabel( "Dimensions: ", gridBox ), 1, 0 );
	QString dims;
	dims.sprintf( "%d x %d x %d pixels", uniVolume->GetDims()[AXIS_X], uniVolume->GetDims()[AXIS_Y], uniVolume->GetDims()[AXIS_Z] );
	gridLayout->addWidget( new QLabel( dims, gridBox ), 1, 1 );

	gridLayout->addWidget( new QLabel( "Total size: ", gridBox ), 2, 0 );
	QString size;
	size.sprintf( "%f x %f x %f mm", uniVolume->Size[AXIS_X], uniVolume->Size[AXIS_Y], uniVolume->Size[AXIS_Z] );
	gridLayout->addWidget( new QLabel( size, gridBox ), 2, 1 );

	gridLayout->addWidget( new QLabel( "Pixel size: ", gridBox ), 3, 0 );
	QString psize;
	psize.sprintf( "%f x %f x %f mm", uniVolume->m_Delta[AXIS_X], uniVolume->m_Delta[AXIS_Y], uniVolume->m_Delta[AXIS_Z] );

	gridLayout->addWidget( new QLabel( psize, gridBox ), 3, 1 );
	} 
      else
	{
	gridLayout->addWidget( new QLabel( "Type: ", gridBox ), 0, 0 );
	gridLayout->addWidget( new QLabel( "Non-uniform rectilinear", gridBox ), 0, 1 );
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
