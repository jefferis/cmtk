/*
//
//  Copyright 2009 SRI International
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

#include <sstream>

void
cmtk::CommandLine::NonOptionParameter
::Evaluate( const size_t argc, const char* argv[], size_t& index )
{
  if ( this->Flag ) 
    *this->Flag = true;
  
  if ( index < argc ) 
    {
    *this->Var = argv[index];
    } 
  else
    {
    throw( Exception( "Argument missing", index ) );
    }
}


mxml_node_t* 
cmtk::CommandLine::NonOptionParameter
::MakeXML( mxml_node_t *const parent, const int index ) const
{
  mxml_node_t *node = NULL;
  if ( this->m_Properties & PROPS_IMAGE )
    {
    node = mxmlNewElement( parent, "image" );
    
    if ( this->m_Properties & PROPS_LABELS )
      mxmlElementSetAttr( node, "type", "label" );
    else
      mxmlElementSetAttr( node, "type", "scalar" );
    }
  else if ( this->m_Properties & PROPS_XFORM )
    {
    node = mxmlNewElement( parent, "transform" );
    mxmlElementSetAttr( node, "fileExtensions", ".txt" );
    }
  else if ( this->m_Properties & PROPS_FILENAME )
    node = mxmlNewElement( parent, "file" );
  else if ( this->m_Properties & PROPS_DIRNAME )
    node = mxmlNewElement( parent, "directory" );
  else 
    node = mxmlNewElement( parent, "string" );

  if ( this->m_Properties & PROPS_OUTPUT )
    mxmlNewText( mxmlNewElement( node, "channel" ), 0, "output" );
  else
    mxmlNewText( mxmlNewElement( node, "channel" ), 0, "input" );

  mxmlNewText( mxmlNewElement( node, "name" ), 0, this->m_Name );
  mxmlNewText( mxmlNewElement( node, "description" ), 0, this->m_Comment );

  if ( index >= 0 )
    {
    std::ostringstream strm;
    strm << index;
    mxmlNewText( mxmlNewElement( node, "index" ), 0, strm.str().c_str() );
    }

  return node;
} 
