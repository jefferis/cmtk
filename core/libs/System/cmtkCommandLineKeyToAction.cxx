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


void
cmtk::CommandLine::KeyToAction
::SetProperties( const long int properties )
{
  this->m_Properties = properties;
}

long int
cmtk::CommandLine::KeyToAction
::GetProperties() const
{
  return this->m_Properties;
}

mxml_node_t*
cmtk::CommandLine::KeyToAction
::MakeXML( mxml_node_t *const node ) const
{
// for some reason Slicer does not accept long options that contain hyphens ("-"), so we replace them.
  std::string xmlKeyStr = this->m_KeyString;
  for ( size_t i = 0; i < xmlKeyStr.length(); ++i )
    {
    if ( xmlKeyStr[i] == '-' )
      xmlKeyStr[i] = '_';
    }

  if ( this->m_Comment.length() )
    {
    mxmlNewText( mxmlNewElement( node, "description" ), 0, this->m_Comment.c_str() );
    }
  
  if ( this->m_KeyString.length() )
    {
    mxmlNewText( mxmlNewElement( node, "name" ), 0, xmlKeyStr.c_str() );
    mxmlNewText( mxmlNewElement( node, "label" ), 0, xmlKeyStr.c_str() );
    }
  if ( this->m_Key )
    {
    const char keyStr[] = { '-', this->m_Key, 0 };
    mxmlNewText( mxmlNewElement( node, "flag" ), 0, keyStr );
    }
  if ( this->m_KeyString.length() )
    {
    mxmlNewText( mxmlNewElement( node, "longflag" ), 0, (std::string( "--" ) + xmlKeyStr).c_str() );
    }

  return node;
}

void
cmtk::CommandLine::KeyToAction
::FormatHelp( std::ostringstream& fmt ) const
{
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

  if ( fmt.str().length() > static_cast<size_t>( CommandLine::HelpTextIndent-2 ) )
    fmt << "\n";
  else
    {
    while ( fmt.str().length() < static_cast<size_t>( CommandLine::HelpTextIndent ) )
      fmt << " ";
    }
  
  if ( this->m_Comment.length() )
    {
    fmt << this->m_Comment;
    }
}

void
cmtk::CommandLine::KeyToAction
::PrintWiki() const
{
  StdOut << "; ";
  if ( this->m_Key )
    {
    StdOut << "<tt>-" << this->m_Key << "</tt>";
    }

  if ( this->m_Key && this->m_KeyString.size() )
    {
    StdOut << " / ";
    }

  if ( this->m_KeyString.size() )
    {
    StdOut << "<tt>--" << this->m_KeyString << "</tt>";
    }

  StdOut << " : ";
  if ( this->m_Comment.length() )
    {
    StdOut << this->m_Comment;
    }
}

bool
cmtk::CommandLine::KeyToAction
::MatchLongOption( const std::string& key ) const
{
  if ( key.length() != this->m_KeyString.length() )
    return false;
  
  for ( size_t i = 0; i < key.length(); ++i )
    {
    if ( (key[i] == '-' || key[i] == '_') && (this->m_KeyString[i] == '-' || this->m_KeyString[i] == '_') )
      continue;

    if ( key[i] != this->m_KeyString[i] )
      return false;
    }
  return true;
}

