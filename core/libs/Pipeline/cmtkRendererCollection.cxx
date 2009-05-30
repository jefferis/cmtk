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

#include <cmtkRendererCollection.h>

#ifdef DEBUG
#  include <stdio.h>
#endif

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

RendererCollection *RendererCollectionInstance = NULL;

RendererCollection::RendererCollection() {};

RendererCollection::~RendererCollection()
{
  while ( !this->m_RendererList.empty() ) 
    {
    this->m_RendererList.back()->Delete();
    this->m_RendererList.pop_back();
    }
}

void RendererCollection::AddRenderer( Renderer *const renderer )
{
#ifdef DEBUG
  if ( ! renderer ) 
    { 
    puts( "Coding error: NULL pointer in AddRenderer" );
    return;
    }
#endif
  renderer->Reference();
  this->m_RendererList.push_back( renderer );
}

void RendererCollection::RemoveRenderer( Renderer *const renderer )
{
#ifdef DEBUG
  if ( ! renderer ) 
    { 
    puts( "Coding error: NULL pointer in RemoveRenderer" );
    return;
    }
#endif
  RendererList::iterator it = this->m_RendererList.begin();
  while ( it != this->m_RendererList.end() ) 
    {
    if ( (*it) == renderer ) 
      {
      (*it)->Delete();
      this->m_RendererList.erase( it );
      }
    ++it;
    }
}

void RendererCollection::Render()
{
  RendererList::iterator it = this->m_RendererList.begin();
  while ( it != this->m_RendererList.end() ) 
    {
    (*it)->Render();
    ++it;
    }
}

int
igsRenderAll()
{
#ifdef DEBUG
  puts( "igsRenderAll() called." ); 
#endif
  if ( RendererCollectionInstance ) 
    {
    RendererCollectionInstance->Render();
    return 1;
    } 
  else
    {
    return 0;
    }
}

} // namespace cmtk
