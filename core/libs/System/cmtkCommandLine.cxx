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

#include <string.h>
#include <sstream>

#include <cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

void 
CommandLine::SetDefaultInfo()
{
  this->m_ProgramInfo[PRG_LCNSE] = "http://www.fsf.org/licensing/licenses/gpl.html";
  this->m_ProgramInfo[PRG_CNTRB] = "Torsten Rohlfing";
  this->m_ProgramInfo[PRG_ACKNL] = "CMTK is supported by the National Institute of Biomedical Imaging and BioEngineering under Grant EB008381";
  this->m_ProgramInfo[PRG_CATEG] = "CMTK.Miscellaneous";
  this->m_ProgramInfo[PRG_DOCUM] = "https://neuro.sri.com/cmtk/wiki/";
  this->m_ProgramInfo[PRG_VERSN] = CMTK_VERSION;
}

CommandLine::KeyActionGroupType::SmartPtr&
CommandLine
::BeginGroup( const char* name, const char* description ) 
{ 
  this->m_KeyActionGroupList.push_back( KeyActionGroupType::SmartPtr( new KeyActionGroupType( name, description ) ) );
  this->m_KeyActionList = &(this->m_KeyActionGroupList.back()->m_KeyActionList);
  return this->m_KeyActionGroupList.back();
}

void
CommandLine
::EndGroup() 
{
  this->m_KeyActionList = &(this->m_KeyActionGroupList.front()->m_KeyActionList);
}

bool
CommandLine::Parse()
{
  this->Index = 1;
  while ( (this->Index < this->ArgC) && (this->ArgV[this->Index][0] == '-') ) 
    {
    // Break at first non-switch argument.
    if ( this->ArgV[this->Index][0] != '-' ) return true;
    
    // Like POSIX, break at "--" terminator.
    if ( !strcmp( this->ArgV[this->Index], "--" ) ) 
      {
      ++this->Index;
      break;
      }
    
    bool found = false;
    if ( this->ArgV[this->Index][1] == '-' ) 
      {
      // long option
      for ( KeyActionListType::iterator it = this->m_KeyActionListComplete.begin(); !found && (it != this->m_KeyActionListComplete.end()); ++it )
	{
	found = (*it)->MatchAndExecute( std::string( this->ArgV[this->Index]+2 ), this->ArgC, this->ArgV, this->Index );
	}
      
      // not found?
      if ( !found ) 
	{
	// Check for "--xml" special option, which produces self description according to Slicer execution model.
	if ( !strcmp( this->ArgV[this->Index], "--xml" ) && !(this->m_Properties & PROPS_NOXML) ) 
	  {
	  this->WriteXML();
	  exit( 0 );
	  }
	
	// Check for "--help" special option, which produces textual description of all command line options
	if ( !strcmp( this->ArgV[this->Index], "--help" ) ) 
	  {
	  this->PrintHelp();
	  exit( 0 );
	  }
	
	// Check for "--wiki" special option, which produces Wiki-markup description of all command line options
	if ( !strcmp( this->ArgV[this->Index], "--wiki" ) ) 
	  {
	  this->PrintWiki();
	  exit( 0 );
	  }
	
	// Check for "--echo" special option, which echoes the command line to stdout. This does not exit the program automatically.
	if ( !strcmp( this->ArgV[this->Index], "--echo" ) ) 
	  {
	  for ( size_t i = 0; i < this->ArgC; ++i )
	    {
	    std::cout << this->ArgV[i] << " ";
	    }
	  std::cout << std::endl;
	  found = true;
	  }
	
	if ( ! found )
	  throw( Exception( std::string("Unknown option: ") + std::string(this->ArgV[this->Index] ) ) );
	}
      } 
    else
      {
      const char* optChar = this->ArgV[this->Index]+1;
      while ( *optChar ) 
	{
	// short option
	for ( KeyActionListType::iterator it = this->m_KeyActionListComplete.begin(); !found && (it != this->m_KeyActionListComplete.end()); ++it )
	  {
	  found = (*it)->MatchAndExecute( *optChar, this->ArgC, this->ArgV, this->Index );
	  }
      
	if ( !found ) 
	  {
	  const char opt[2] = { *optChar, 0 };
	  throw( Exception( std::string("Unknown option: -") + std::string(opt) ) );
	  }
	
	++optChar; // next short option in this block, POSIX style.
	}
      }
    
    ++this->Index;
    } // while

  for ( NonOptionParameterListType::iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it, ++this->Index )
    {
    if ( this->Index >= this->ArgC )
      {     
      if ( ! ((*it)->m_Properties & PROPS_OPTIONAL) )
	throw( Exception( "Insufficient number of command line arguments", this->Index ) );
      }
    else
      {
      (*it)->Evaluate( this->ArgC, this->ArgV, this->Index );
      }
    }
  
  for ( NonOptionParameterVectorListType::iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it, ++this->Index )
    {
    if ( this->Index >= this->ArgC )
      {     
      if ( ! ((*it)->m_Properties & PROPS_OPTIONAL) )
	throw( Exception( "Insufficient number of command line arguments", this->Index ) );
      }
    else
      {
      (*it)->Evaluate( this->ArgC, this->ArgV, this->Index );
      }
    }
  
  return true;
}

void
CommandLine::PrintHelp
() const
{
  const size_t lineWidth = StdErr.GetLineWidth();

  ProgramPropertiesMapType::const_iterator ppit = this->m_ProgramInfo.find(PRG_TITLE);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdErr << "TITLE:\n\n";
    StdErr.FormatText( ppit->second, 5 ) << "\n";
    }

  ppit = this->m_ProgramInfo.find(PRG_DESCR);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdErr << "\nDESCRIPTION:\n\n";
    StdErr.FormatText( ppit->second, 5 ) << "\n";
    }

  ppit = this->m_ProgramInfo.find(PRG_SYNTX);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdErr << "\nSYNTAX:\n\n";
    StdErr.FormatText( ppit->second, 5 ) << "\n";
    }
  else
    {
    if ( this->m_NonOptionParameterList.size() || this->m_NonOptionParameterVectorList.size() )
      {
      StdErr << "\nSYNTAX:\n\n";

      std::ostringstream fmt;
      fmt << "[options] ";
      for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	{
	fmt << (*it)->m_Name << " ";
	}
      for ( NonOptionParameterVectorListType::const_iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it )
	{
	fmt << (*it)->m_Name << " ";
	}
      StdErr.FormatText( fmt.str(), 5, lineWidth );

      StdErr << "\n  where\n";

      const int indent = 20;
      for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	{
	fmt.str("");

	StdErr << "\n";
	fmt << (*it)->m_Name << " = ";
	if ( fmt.str().length() > static_cast<size_t>( indent-2 ) )
	  fmt << "\n";
	else
	  {
	  while ( fmt.str().length() < static_cast<size_t>( indent ) )
	    fmt << " ";
	  }
	fmt << (*it)->m_Comment;
	StdErr.FormatText( fmt.str(), 5+indent, lineWidth, -indent ) << "\n";
	}

      for ( NonOptionParameterVectorListType::const_iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it )
	{
	fmt.str("");
	
	StdErr << "\n";
	fmt << (*it)->m_Name << " = ";
	if ( fmt.str().length() > static_cast<size_t>( indent-2 ) )
	  fmt << "\n";
	else
	  {
	  while ( fmt.str().length() < static_cast<size_t>( indent ) )
	    fmt << " ";
	  }
	fmt << (*it)->m_Comment;
	StdErr.FormatText( fmt.str(), 5+indent, lineWidth, -indent ) << "\n";
	}
      }
    }

  StdErr << "\nLIST OF SUPPORTED OPTIONS:\n\n";

  for ( KeyActionGroupListType::const_iterator grp = this->m_KeyActionGroupList.begin(); grp != this->m_KeyActionGroupList.end(); ++grp )
    {
    const std::string& name = (*grp)->m_Name;

    size_t indent = 0;
    if ( name != "MAIN" )
      {
      StdErr << (*grp)->m_Description << "\n\n";
      indent = 2;
      }

    const KeyActionListType& kal = (*grp)->m_KeyActionList;
    for ( KeyActionListType::const_iterator it = kal.begin(); it != kal.end(); ++it )
      {
      (*it)->PrintHelp( indent );
      }
    }
  
  StdErr << "\n";
}

void
CommandLine::PrintWiki
() const
{
  ProgramPropertiesMapType::const_iterator ppit = this->m_ProgramInfo.find(PRG_TITLE);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << "== Title ==\n\n";
    StdOut << ppit->second << "\n\n";
    }

  ppit = this->m_ProgramInfo.find(PRG_DESCR);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << "== Description ==\n\n";
    StdOut << ppit->second << "\n\n";
    }
  
  ppit = this->m_ProgramInfo.find(PRG_SYNTX);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << "== Syntax ==\n\n";
    StdOut << ppit->second << "\n\n";
    }
  else
    {
    if ( this->m_NonOptionParameterList.size() || this->m_NonOptionParameterVectorList.size() )
      {
      StdOut << "== Syntax ==\n\n";
      
      StdOut << ": <tt>[options] ";
      for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	{
	StdOut << (*it)->m_Name << " ";
	}
      for ( NonOptionParameterVectorListType::const_iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it )
	{
	StdOut << (*it)->m_Name << " ";
	}
      StdOut << "</tt>\n\nwhere\n";

      for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	{
	StdOut << "\n";
	StdOut << "; <tt>" << (*it)->m_Name << "</tt> : ";
	StdOut << (*it)->m_Comment << "\n";;
	}
      for ( NonOptionParameterVectorListType::const_iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it )
	{
	StdOut << "\n";
	StdOut << "; <tt>" << (*it)->m_Name << "</tt> : ";
	StdOut << (*it)->m_Comment << "\n";;
	}
      }
    }
  
  StdOut << "\n== List of Supported Options ==\n\n";

  for ( KeyActionGroupListType::const_iterator grp = this->m_KeyActionGroupList.begin(); grp != this->m_KeyActionGroupList.end(); ++grp )
    {
    const std::string& name = (*grp)->m_Name;
    if ( name != "MAIN" )
      {
      StdOut << "=== " << (*grp)->m_Description << " ===\n\n";
      }

    const KeyActionListType& kal = (*grp)->m_KeyActionList;
    for ( KeyActionListType::const_iterator it = kal.begin(); it != kal.end(); ++it )
      {
      (*it)->PrintWikiWithPrefix();
      StdOut << "\n";
      }
    }
  
  StdOut << "\n";
}

Console& operator<<( Console& console, CommandLine::Exception e )
{
  console << e.Message << " [argument #" << e.Index << "]\n";
  return console;
}

} // namespace cmtk
