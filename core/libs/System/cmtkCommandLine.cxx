/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include "cmtkCommandLine.h"

#include <System/cmtkConsole.h>
#include <System/cmtkExitException.h>
#include <System/cmtkDebugOutput.h>
#include <System/cmtkThreads.h>

#include <string.h>
#include <sstream>

namespace
cmtk
{

/** \addtogroup System */
//@{

CommandLine::CommandLine( const int properties ) 
  : ArgC( 0 ),
    ArgV( NULL ),
    m_Properties( properties )
{
  this->SetDefaultInfo();    
  this->BeginGroup( "MAIN", "Main Options" );
}

CommandLine::~CommandLine()
{
  if ( this->Index < this->ArgC-1 )
    {
    StdErr << "WARNING: the following command line arguments were not used: \n";
    for ( size_t i = this->Index; i < this->ArgC; ++i )
      {
      StdErr << this->ArgV[i] << " ";
      }
    StdErr << "\n";
    }
}

void 
CommandLine::SetDefaultInfo()
{
  this->m_ProgramInfo[PRG_LCNSE] = "http://www.fsf.org/licensing/licenses/gpl.html";
  this->m_ProgramInfo[PRG_CNTRB] = "Torsten Rohlfing, Michael P. Hasak, Greg Jefferis, Calvin R. Maurer, Daniel B. Russakoff";
  this->m_ProgramInfo[PRG_ACKNL] = "CMTK is supported by the National Institute of Biomedical Imaging and BioEngineering under Grant EB008381";
  this->m_ProgramInfo[PRG_CATEG] = "CMTK.Miscellaneous";
  this->m_ProgramInfo[PRG_DOCUM] = "https://neuro.sri.com/cmtk/wiki/";
  this->m_ProgramInfo[PRG_VERSN] = CMTK_VERSION_STRING;

  this->BeginGroup( "GLOBAL", "Global Toolkit Options" )->SetProperties( Self::PROPS_NOXML );
  this->AddCallback( Self::Key( "help" ), &Self::CallbackInternal, "Write list of command line options to standard output." );
  this->AddCallback( Self::Key( "wiki" ), &Self::CallbackInternal, "Write list of command line options to standard output in MediaWiki markup." );

  if (! (this->m_Properties & PROPS_NOXML) ) 
    this->AddCallback( Self::Key( "xml" ), &Self::CallbackInternal, "Write command line syntax specification in XML markup (for Slicer integration)." );

  this->AddCallback( Self::Key( "version" ), &Self::CallbackInternal, "Write toolkit version to standard output." );
  this->AddCallback( Self::Key( "echo" ), &Self::CallbackInternal, "Write the current command line to standard output." );
  this->AddCallback( Self::Key( "verbose-level" ), &DebugOutput::SetGlobalLevel, "Set verbosity level." );
  this->AddCallback( Self::Key( 'v', "verbose" ), &DebugOutput::IncGlobalLevel, "Increment verbosity level by 1 (deprecated; supported for backward compatibility)." );
  this->AddCallback( Self::Key( "threads" ), &Threads::SetNumberOfThreads, "Set maximum number of parallel threads (for POSIX threads and OpenMP)." );
  this->EndGroup();
}

void 
CommandLine::CallbackInternal()
{
  StdErr << "ERROR: cmtk::CommandLine::CallbackInternal should never be called.\n";
  throw ExitException( 1 );
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
CommandLine::Parse( const int argc, const char* argv[] ) throw( ExitException, Self::Exception )
{
  this->ArgC = argc;
  this->ArgV = argv;

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
      // Check for "--version" special option, which prints the CMTK version that this tool is part of.
      if ( !strcmp( this->ArgV[this->Index], "--version" ) ) 
	{
	StdOut << this->m_ProgramInfo[PRG_VERSN] << "\n";
	throw ExitException( 0 );
	}
      
      // Check for "--xml" special option, which produces self description according to Slicer execution model.
      if ( !strcmp( this->ArgV[this->Index], "--xml" ) && !(this->m_Properties & PROPS_NOXML) ) 
	{
	this->WriteXML();
	throw ExitException( 0 );
	}
      
      // Check for "--help" special option, which produces textual description of all command line options
      if ( !strcmp( this->ArgV[this->Index], "--help" ) ) 
	{
	this->PrintHelp();
	throw ExitException( 0 );
	}
      
      // Check for "--wiki" special option, which produces Wiki-markup description of all command line options
      if ( !strcmp( this->ArgV[this->Index], "--wiki" ) ) 
	{
	this->PrintWiki();
	throw ExitException( 0 );
	}
      
      // Check for "--echo" special option, which echoes the command line to stderr. This does not exit the program automatically.
      if ( !strcmp( this->ArgV[this->Index], "--echo" ) ) 
	{
	for ( size_t i = 0; i < this->ArgC; ++i )
	  {
	  std::cerr << this->ArgV[i] << " ";
	  }
	std::cerr << std::endl;
	found = true;
	}
	
      if ( !found )
	{
	for ( KeyActionListType::iterator it = this->m_KeyActionListComplete.begin(); !found && (it != this->m_KeyActionListComplete.end()); ++it )
	  {
	  found = (*it)->MatchAndExecute( std::string( this->ArgV[this->Index]+2 ), this->ArgC, this->ArgV, this->Index );
	  }
	}
      
      if ( ! found )
	throw( Exception( std::string("Unknown option: ") + std::string( this->ArgV[this->Index] ) ) );
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
  const size_t lineWidth = StdOut.GetLineWidth();

  ProgramPropertiesMapType::const_iterator ppit = this->m_ProgramInfo.find(PRG_TITLE);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << "TITLE:\n\n";
    StdOut.FormatText( ppit->second, 5 ) << "\n";
    }

  ppit = this->m_ProgramInfo.find(PRG_DESCR);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << "\nDESCRIPTION:\n\n";
    StdOut.FormatText( ppit->second, 5 ) << "\n";
    }

  ppit = this->m_ProgramInfo.find(PRG_SYNTX);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << "\nSYNTAX:\n\n";
    StdOut.FormatText( ppit->second, 5 ) << "\n";
    }
  else
    {
    if ( this->m_NonOptionParameterList.size() || this->m_NonOptionParameterVectorList.size() )
      {
      StdOut << "\nSYNTAX:\n\n";

      std::ostringstream fmt;
      fmt << "[options] ";
      for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	{
	if ( ! (*it)->m_Comment.empty() )
	  {
	  fmt << (*it)->m_Name << " ";
	  }
	}
      for ( NonOptionParameterVectorListType::const_iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it )
	{
	if ( ! (*it)->m_Comment.empty() )
	  {
	  fmt << (*it)->m_Name << " ";
	  }
	}
      StdOut.FormatText( fmt.str(), 5, lineWidth );

      StdOut << "\n  where\n";

      const int indent = 20;
      for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	{
	if ( ! (*it)->m_Comment.empty() )
	  {
	  fmt.str("");
	  
	  StdOut << "\n";
	  fmt << (*it)->m_Name << " = ";
	  if ( fmt.str().length() > static_cast<size_t>( indent-2 ) )
	    fmt << "\n";
	  else
	    {
	    while ( fmt.str().length() < static_cast<size_t>( indent ) )
	      fmt << " ";
	    }
	  fmt << (*it)->m_Comment;
	  StdOut.FormatText( fmt.str(), 5+indent, lineWidth, -indent ) << "\n";
	  }
	}

      for ( NonOptionParameterVectorListType::const_iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it )
	{
	if ( ! (*it)->m_Comment.empty() )
	  {
	  fmt.str("");
	  
	  StdOut << "\n";
	  fmt << (*it)->m_Name << " = ";
	  if ( fmt.str().length() > static_cast<size_t>( indent-2 ) )
	    fmt << "\n";
	  else
	    {
	    while ( fmt.str().length() < static_cast<size_t>( indent ) )
	      fmt << " ";
	    }
	  fmt << (*it)->m_Comment;
	  StdOut.FormatText( fmt.str(), 5+indent, lineWidth, -indent ) << "\n";
	  }
	}
      }
    }

  StdOut << "\nLIST OF SUPPORTED OPTIONS:\n\n";

  for ( KeyActionGroupListType::const_iterator grp = this->m_KeyActionGroupList.begin(); grp != this->m_KeyActionGroupList.end(); ++grp )
    {
    if ( ! (*grp)->m_KeyActionList.empty() )
      {
      StdOut << (*grp)->m_Description << "\n\n";
      
      const size_t indent = 2;
      
      const KeyActionListType& kal = (*grp)->m_KeyActionList;
      for ( KeyActionListType::const_iterator it = kal.begin(); it != kal.end(); ++it )
	{
	(*it)->PrintHelp( indent );
	}
      }
    }
  
  StdOut << "\n";
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
    if ( ! (*grp)->m_KeyActionList.empty() )
      {
      StdOut << "=== " << (*grp)->m_Description << " ===\n\n";
      
      const KeyActionListType& kal = (*grp)->m_KeyActionList;
      for ( KeyActionListType::const_iterator it = kal.begin(); it != kal.end(); ++it )
	{
	(*it)->PrintWikiWithPrefix();
	StdOut << "\n";
	}
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
