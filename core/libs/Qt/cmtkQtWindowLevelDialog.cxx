/*
//
//  Copyright 1997-2009 Torsten Rohlfing
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

#include "cmtkQtWindowLevelDialog.h"

#include "Qt/cmtkQtIcons.h"

#include <qlayout.h>
#include <qboxlayout.h>

namespace
cmtk
{

/** \addtogroup Qt */
//@{

QtWindowLevelDialog::QtWindowLevelDialog
( QWidget* parent, bool modal, Qt::WFlags f )
  : QDialog( parent, f )
{
  this->setModal( modal );
  this->setWindowIcon( QtIcons::WindowIcon() );
  this->setWindowTitle( "Window/Level Control" );

  QBoxLayout* layout = new QVBoxLayout( this );

  Controls = new QtWindowLevelControls( this );
  QObject::connect( Controls, SIGNAL( colormap( Study::SmartPtr& ) ), SIGNAL( colormapChanged( Study::SmartPtr& ) ) );
  layout->addWidget( Controls );
}

void
QtWindowLevelDialog::slotSetStudy( Study::SmartPtr& study )
{
  Controls->slotSetStudy( study );
}

} // namespace cmtk
