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
cmtk::CommandLine::KeyToActionSingle
::MakeXML( mxml_node_t *const parent ) const
{
  return this->Superclass::MakeXML( this->m_Action->MakeXML( parent ) );
}

void
cmtk::CommandLine::KeyToActionSingle
::PrintHelp( const size_t globalIndent ) const
{
  std::ostringstream fmt;
  this->Superclass::FormatHelp( fmt );

  this->m_Action->PrintHelp( fmt );
  StdErr.FormatText( fmt.str(), CommandLine::HelpTextIndent + globalIndent, 80, -CommandLine::HelpTextIndent ) << "\n";  
}

bool
cmtk::CommandLine::KeyToActionSingle
::MatchAndExecute( const std::string& key, const size_t argc, const char* argv[], size_t& index )
{
  if ( this->MatchLongOption( std::string( key ) ) )
    {
    this->m_Action->Evaluate( argc, argv, index );
    return true;
    }
  return false;
}

bool
cmtk::CommandLine::KeyToActionSingle
::MatchAndExecute( const char keyChar, const size_t argc, const char* argv[], size_t& index )
{
  if ( this->m_Key == keyChar )
    {
    this->m_Action->Evaluate( argc, argv, index );
    return true;
    }

  return false;
}
