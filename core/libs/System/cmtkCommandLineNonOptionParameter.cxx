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
    if ( ! (this->m_Properties & PROPS_OPTIONAL) )
      throw( Exception( "Argument missing", index ) );
    }
}


mxml_node_t* 
cmtk::CommandLine::NonOptionParameter
::MakeXML( mxml_node_t *const parent, const int index ) const
{
  mxml_node_t *node = Item::Helper<const char*>::MakeXML( this, parent );

  if ( node )
    {
    if ( this->m_Name )
      {
      mxmlNewText( mxmlNewElement( node, "name" ), 0, this->m_Name );
      mxmlNewText( mxmlNewElement( node, "label" ), 0, this->m_Name );
      }
    
    if ( this->m_Comment )
      {
      mxmlNewText( mxmlNewElement( node, "description" ), 0, this->m_Comment );
      }
    
    if ( index >= 0 )
      {
      std::ostringstream strm;
      strm << index;
      mxmlNewText( mxmlNewElement( node, "index" ), 0, strm.str().c_str() );
      }
    }

  return node;
} 

std::string
cmtk::CommandLine::NonOptionParameter
::GetParamTypeString() const
{
  return Item::Helper<const char*>::GetParamTypeString( this );
}
