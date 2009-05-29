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
#include <cmtkCommandLine.h>

mxml_node_t* 
cmtk::CommandLine::Callback
::MakeXML(  mxml_node_t *const parent ) const 
{
  mxml_node_t* node = NULL;
  if ( this->m_FuncArg )
    {
    node = mxmlNewElement( parent, "boolean" );
    mxml_node_t *dflt = mxmlNewElement( node, "default" );
    mxmlNewText( dflt, 0, "false" );
    }
  else if ( this->m_FuncArg )
    {
    node = mxmlNewElement( parent, "string" );
    }
  else if ( this->m_FuncIntArg )
    {
    node = mxmlNewElement( parent, "integer" );
    }
  else if ( this->m_FuncDblArg )
    {
    node = mxmlNewElement( parent, "double" );
    }
  else if ( this->m_FuncMultiArg )
    {
    node = mxmlNewElement( parent, "string-vector" );
    }
  return node;
}
