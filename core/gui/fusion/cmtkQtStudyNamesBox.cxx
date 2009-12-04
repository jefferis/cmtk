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

#include <cmtkQtStudyNamesBox.h>
#include <QHBoxLayout>
#include <QLabel>

namespace
cmtk
{

QtStudyNamesBox::QtStudyNamesBox( QWidget *parent, const char* name )
  : QWidget( parent, name )
{
  Layout = new QHBoxLayout( this );

  Label = new QLabel( this );
  Layout->addWidget( Label );

  Layout->addItem( new QSpacerItem( 10, 0, QSizePolicy::Fixed ) );

  ComboBox = new QComboBox( this );
  Layout->addWidget( ComboBox );
  QObject::connect( ComboBox, SIGNAL( activated( int ) ), SIGNAL( signalActivated( int ) ) );
  QObject::connect( ComboBox, SIGNAL( activated( const QString& ) ), SIGNAL( signalActivated( const QString& ) ) );

  Layout->addItem( new QSpacerItem( 1, 1 ) );
}

void
QtStudyNamesBox::slotUpdateStudySelection( const QStringList *namesList )
{
  // find out where (if) the current left study is in the new list so we can
  // switch to it after the change
  int switchToStudy = 0;
  if ( ComboBox->count() ) 
    {
    QString studyName = ComboBox->currentText();
    switchToStudy = namesList->findIndex( studyName );
    }
  ComboBox->clear();
  ComboBox->insertStringList( *namesList );
  ComboBox->setCurrentItem( switchToStudy );
}

void
QtStudyNamesBox::slotSetCurrentText( const QString& text )
{
  ComboBox->setCurrentText( text );
}

} // namespace cmtk
