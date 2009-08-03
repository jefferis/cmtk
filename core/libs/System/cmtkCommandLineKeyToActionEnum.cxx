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

void
cmtk::CommandLine::KeyToActionEnum
::PrintHelp( const size_t globalIndent ) const
{
  std::ostringstream fmt;
  this->Superclass::FormatHelp( fmt );

  fmt << "\nSupported values: ";
  for ( EnumGroupBase::const_iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
    {
    fmt << "\"" << (*it)->m_KeyString << "\", ";
    }
  
  const std::string defaultKey = this->m_EnumGroup->GetDefaultKey();
  if ( defaultKey.length() )
    {
    fmt << "where the default is \"" << defaultKey << "\", ";
    }
  
  fmt << "or use one of the following";
  
  StdErr.FormatText( fmt.str(), CommandLine::HelpTextIndent + globalIndent, 80, -CommandLine::HelpTextIndent ) << "\n";  
  
  for ( EnumGroupBase::const_iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
    {
    (*it)->PrintHelp( globalIndent + 10 );
    }
}

void
cmtk::CommandLine::KeyToActionEnum
::PrintWiki() const
{
  std::ostringstream fmt;
  this->Superclass::FormatHelp( fmt );

  StdOut << "\nSupported values: ";
  for ( EnumGroupBase::const_iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
    {
    StdOut << "\"" << (*it)->m_KeyString << "\", ";
    }
  
  const std::string defaultKey = this->m_EnumGroup->GetDefaultKey();
  if ( defaultKey.length() )
    {
    StdOut << "where the default is \"" << defaultKey << "\", ";
    }
  
  StdOut << "or use one of the following\n";
  
  for ( EnumGroupBase::const_iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
    {
    (*it)->PrintWiki();
    }
}

mxml_node_t*
cmtk::CommandLine::KeyToActionEnum
::MakeXML( mxml_node_t *const parent ) const
{
  mxml_node_t *node = mxmlNewElement( parent, "string-enumeration" );
  
  mxml_node_t* defaultElement = mxmlNewElement( node, "default" );
  mxmlNewText( defaultElement, 0, this->m_EnumGroup->GetDefaultKey().c_str() );
  
  for ( EnumGroupBase::const_iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
    {      
    mxml_node_t* element = mxmlNewElement( node, "element" );
    mxmlNewText( element, 0, (*it)->m_KeyString.c_str() );
    }
  
  return this->Superclass::MakeXML( node );
}

bool
cmtk::CommandLine::KeyToActionEnum
::MatchAndExecute( const std::string& key, const size_t argc, const char* argv[], size_t& index )
{
  if ( this->MatchLongOption( std::string( key ) ) )
    {
    // check if optional argument matches suboption
    if ( this->m_EnumGroup )
      {      
      for ( EnumGroupBase::iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
	{
	size_t ii = index+1;
	if ( (*it)->MatchAndExecute( argv[ii], argc, argv, ii ) )
	  {
	  index = ii;
	  return true;
	  }
	}
      }
    }

  // also check for direct matches with the suboptions
  if ( this->m_EnumGroup )
    {
    for ( EnumGroupBase::iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
      {
      if ( (*it)->MatchAndExecute( key, argc, argv, index ) )
	{
	return true;
	}
      }
    }
  
  return false;
}

bool
cmtk::CommandLine::KeyToActionEnum
::MatchAndExecute( const char keyChar, const size_t argc, const char* argv[], size_t& index )
{
  // check if optional argument matches suboption
  for ( EnumGroupBase::iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
    {
    size_t ii = index+1;
    if ( (*it)->MatchAndExecute( argv[ii], argc, argv, ii ) )
      {
      index = ii;
      return true;
      }
    }
  
  // also check for direct matches with the suboptions
  for ( EnumGroupBase::iterator it = this->m_EnumGroup->begin(); it != this->m_EnumGroup->end(); ++it )
    {
    if ( (*it)->MatchAndExecute( keyChar, argc, argv, index ) )
      {
      return true;
      }
    }

  return false;
}
