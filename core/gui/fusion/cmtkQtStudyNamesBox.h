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

#ifndef __cmtkQtStudyNamesBox_h_included_
#define __cmtkQtStudyNamesBox_h_included_

#include <cmtkconfig.h>

#include <qwidget.h>
#include <qlayout.h>
#include <qlabel.h>
#include <qstring.h>
#include <qstringlist.h>
#include <qcombobox.h>
#include <Q3HBoxLayout>

namespace
cmtk
{

/// Combo box for study names.
class QtStudyNamesBox :
  /// This is a widget.
  public QWidget 
{
  Q_OBJECT

public:
  /// Constructor.
  QtStudyNamesBox( QWidget *parent, const char *name );

  /// Get currently selected name.
  const QString GetCurrentName() const {
    return ComboBox->currentText();
  }

signals:
  /// Forwarded signal from ComboBox.
  void signalActivated( int );

  /// Forwarded signal from ComboBox.
  void signalActivated( const QString& );

public slots:
  /// Set box label.
  void slotSetLabel( const QString& label ) {
    Label->setText( label );
  }

  /// Explicitly set selected text.
  void slotSetCurrentText( const QString& text );

  /// Update list of study names.
  void slotUpdateStudySelection( const QStringList *namesList );

private:
  /// The enclosing layout.
  Q3HBoxLayout* Layout;

  /// The combo box widget.
  QComboBox* ComboBox;

  /// The label.
  QLabel* Label;
};

} // namespace cmtk

#endif // #ifndef __cmtkQtStudyNamesBox_h_included_
