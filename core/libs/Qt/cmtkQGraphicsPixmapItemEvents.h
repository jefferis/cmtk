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

#ifndef __cmtkQGraphicsPixmapItemEvents_h_included_
#define __cmtkQGraphicsPixmapItemEvents_h_included_

#include <cmtkconfig.h>

#include <QtCore/QObject>

#include <QtGui/QGraphicsPixmapItem>
#include <QtGui/QGraphicsSceneMouseEvent>

namespace
cmtk
{

/// Class that derives from Qt's pixmap graphics item and signals events it receives.
class QGraphicsPixmapItemEvents : public QObject, public QGraphicsPixmapItem
{
  Q_OBJECT
  
public:
  /// This class.
  typedef QGraphicsPixmapItemEvents Self;

  /// Parent class.
  typedef QGraphicsPixmapItem Superclass;

  /// Default constructor.
  QGraphicsPixmapItemEvents( QGraphicsItem* parent = 0 ) : QGraphicsPixmapItem ( parent ) {}

  /// Constructor.
  QGraphicsPixmapItemEvents( const QPixmap& pixmap, QGraphicsItem* parent = 0 ) : QGraphicsPixmapItem ( pixmap, parent ) {}

signals:
  /// Signal that is sent when "mouse press" event is received.
  void mousePressed( QGraphicsSceneMouseEvent* event );

protected:
  /// Catch mouse press event.
  virtual void mousePressEvent( QGraphicsSceneMouseEvent* event )
  {
    emit( mousePressed( event ) );
    Superclass::mousePressEvent( event );
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkQGraphicsPixmapItemEvents_h_included_
