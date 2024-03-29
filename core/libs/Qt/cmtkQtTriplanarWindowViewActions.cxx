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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include "cmtkQtTriplanarWindow.h"

namespace
cmtk
{

void 
QtTriplanarWindow::slotView25()
{
  this->slotSetZoom( 25 );
}

void 
QtTriplanarWindow::slotView33()
{
  this->slotSetZoom( 33 );
}

void 
QtTriplanarWindow::slotView50()
{
  this->slotSetZoom( 50 );
}

void 
QtTriplanarWindow::slotView100()
{
  this->slotSetZoom( 100 );
}

void 
QtTriplanarWindow::slotView200()
{
  this->slotSetZoom( 200 );
}

void 
QtTriplanarWindow::slotView300()
{
  this->slotSetZoom( 300 );
}

void 
QtTriplanarWindow::slotView400()
{
  this->slotSetZoom( 400 );
}

void 
QtTriplanarWindow::slotView500()
{
  this->slotSetZoom( 500 );
}

void 
QtTriplanarWindow::slotViewInterpolation()
{
  this->slotRenderAll();
}

void
QtTriplanarWindow::slotViewCrosshair()
{
  const bool mode = this->m_CrosshairAction->isChecked();
  ScrollRenderViewAx->GetRenderImage()->SetCrosshairMode( mode );
  ScrollRenderViewCo->GetRenderImage()->SetCrosshairMode( mode );
  ScrollRenderViewSa->GetRenderImage()->SetCrosshairMode( mode );
  this->slotRenderAll();
}
      
void
QtTriplanarWindow::slotViewCheckerbox()
{
  this->slotRenderAll();
}
      
} // namespace cmtk
