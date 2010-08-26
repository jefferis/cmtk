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

#include "Pipeline/cmtkRenderer.h"

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

Renderer::Renderer()
{
  Input = NULL;
  Active = true;
  RenderPending = false;
}

Renderer::~Renderer()
{
  if ( Input != NULL ) Input->Delete();
}

void
Renderer::SetInput( ImageRGB *const input )
{
  ReplaceObject( Input, input );
};

long
Renderer::Update()
{
  if ( this->IsActive() )
    this->CheckInputForUpdate( Input );
  return this->Superclass::Update();
}

void 
Renderer::Render()
{
  // Is this renderer being updated already, ie. do we have a recursion here?
  // If no: go through it.
  if ( ! RenderPending ) {
    RenderPending = true;

    // Fake modification to make sure we WILL update ourselves.
    this->UpdateModifiedTime();
    this->Update();

    RenderPending = false;
  }
}

} // namespace cmtk
