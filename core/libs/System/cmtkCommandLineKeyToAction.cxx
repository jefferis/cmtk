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

#include <sstream>

mxml_node_t*
cmtk::CommandLine::KeyToAction::MakeXML( mxml_node_t *const parent ) const
{
  mxml_node_t *node = this->m_Action->MakeXML( parent );
  if ( this->m_Comment )
    {
    mxmlNewText( mxmlNewElement( node, "description" ), 0, this->m_Comment );
    }
  if ( this->m_Key )
    {
    const char keyStr[] = { '-', this->m_Key, 0 };
    mxmlNewText( mxmlNewElement( node, "flag" ), 0, keyStr );
    }
  if ( this->m_KeyString.size() )
    {	
    mxmlNewText( mxmlNewElement( node, "longflag" ), 0, (std::string( "--" ) + this->m_KeyString).c_str() );
    }

  return node;
}

void
cmtk::CommandLine::KeyToAction::PrintHelp( const size_t globalIndent ) const
{
  std::ostringstream fmt;
  if ( this->m_Key )
    {
    fmt << "-" << this->m_Key;
    }

  if ( this->m_Key && this->m_KeyString.size() )
    {
    fmt << " / ";
    }

  if ( this->m_KeyString.size() )
    {
    fmt << "--" << this->m_KeyString;
    }

  const int indent = 20;
  if ( fmt.str().length() > static_cast<size_t>( indent-2 ) )
    fmt << "\n";
  else
    {
    while ( fmt.str().length() < static_cast<size_t>( indent ) )
      fmt << " ";
    }
  
  if ( this->m_Comment )
    {
    fmt << this->m_Comment;
    }
  
  StdErr.FormatText( fmt.str(), indent + globalIndent, 80, -indent ) << "\n";
}
