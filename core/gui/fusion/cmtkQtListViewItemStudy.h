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

#ifndef __cmtkQtListViewItemStudy_h_included_
#define __cmtkQtListViewItemStudy_h_included_

#include <cmtkconfig.h>

#include <q3listview.h>

#include <cmtkMacros.h>
#include <cmtkStudy.h>

namespace
cmtk
{

/** Class that extents QListViewItem with an Study payload.
 */
class QtListViewItemStudy :
  /// Inherit from Qt's list view item
  public Q3ListViewItem
{
  /// Pointer to the Study object.
  cmtkGetSetMacro(Study::SmartPtr,Study);

public:
  /// Constructor when adding as root item.
  QtListViewItemStudy( Q3ListView *const parent, Study::SmartPtr& study ) 
    : Q3ListViewItem( parent, study->GetName(), study->GetDescription() ),
      m_Study(NULL)
  {
    this->SetStudy( study );
  }

  /// Constructor when adding below another list view item.
  QtListViewItemStudy
  ( Q3ListViewItem *const parent, Study::SmartPtr& study )
    : Q3ListViewItem( parent, study->GetName(), study->GetDescription() ),
      m_Study(NULL)
  {
    this->SetStudy( study );
  }

  virtual ~QtListViewItemStudy() 
  {
    this->SetStudy( Study::SmartPtr::Null );
  }
};

} // namespace cmtk

#endif // #ifndef __cmtkQtListViewItemStudy_h_included_
