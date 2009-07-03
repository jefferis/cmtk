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

void 
CommandLine::Callback::Evaluate
( const size_t argc, const char* argv[], size_t& index )
{
  // callback with argument?
  if ( this->m_FuncArg ) 
    {
    if ( index+1 < argc ) 
      {
      const char* error = this->m_FuncArg( argv[index+1] );
      if ( error ) 
	{
	throw( Exception( error, index ) );
	} 
      else 
	{
	++index;
	}
      } 
    else
      {
      throw( Exception( "Option needs an argument.", index ) );
      }
    } 
  else
    // callback with integer argument?
    if ( this->m_FuncIntArg ) 
      {
      if ( index+1 < argc ) 
	{
	const char* error = this->m_FuncIntArg( ConvertStrToLong( argv[index+1] ) );
	if ( error ) 
	  {
	  throw( Exception( error, index ) );
	  } 
	else 
	  {
	  ++index;
	  }
	} 
      else
	{
	throw( Exception( "Option needs an integer argument.", index ) );
	}
      } 
  else
    // callback with double argument?
    if ( this->m_FuncDblArg ) 
      {
      if ( index+1 < argc ) 
	{
	const char* error = this->m_FuncDblArg( ConvertStrToDouble( argv[index+1] ) );
	if ( error ) 
	  {
	  throw( Exception( error, index ) );
	  } 
	else 
	  {
	  ++index;
	  }
	} 
      else
	{
	throw( Exception( "Option needs a floating point argument.", index ) );
	}
      } 
  else
    // multiple arguments to callback?
    if ( this->m_FuncMultiArg ) 
      {
      if ( index+1 < argc ) 
	{
	int argsUsed = 0;
	const char* error = this->m_FuncMultiArg( argv+index+1, argsUsed );
	if ( error ) 
	  {
	  throw( Exception( error, index ) );
	  } 
	else 
	  {
	  index += argsUsed;
	  }
	} 
      else
	{
	throw( Exception( "Option needs an argument", index ) );
	}
      } 
    else
      {
      // no argument to callback
      const char* error = this->m_Func();
      if ( error ) 
	{
	throw( Exception( error, index ) );
	}
      }
}

bool
CommandLine::MatchLongOption( const std::string& s1, const std::string& s2 ) const
{
  if ( s1.length() != s2.length() )
    return false;

  for ( size_t i = 0; i < s1.length(); ++i )
    {
    if ( (s1[i] == '-' || s1[i] == '_') && (s2[i] == '-' || s2[i] == '_') )
      continue;

    if ( s1[i] != s2[i] )
      return false;
    }
  return true;
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
    
    SmartPointer<Item> item(NULL);
    if ( this->ArgV[this->Index][1] == '-' ) 
      {
      // long option
      for ( KeyActionListType::iterator it = this->m_KeyActionListComplete.begin(); it != this->m_KeyActionListComplete.end(); ++it )
	{
	if ( MatchLongOption( (*it)->m_KeyString, std::string( this->ArgV[this->Index]+2 ) ) )
	  {
	  item = (*it)->m_Action;
	  }
	}
      
      // not found?
      if ( !item ) 
	{
	// Check for "--xml" special option, which produces self description according to Slicer execution model.
	if ( !strcmp( this->ArgV[this->Index], "--xml" ) ) 
	  {
	  this->WriteXML();
	  exit( 0 );
	  }
	
	// Check for "--help" special option, which produces textual description of all command line options
	if ( !strcmp( this->ArgV[this->Index], "--help" ) ) 
	  {
	  this->PrintHelp();
	  exit( 2 );
	  }
	
	throw( Exception( std::string("Unknown option: ") + std::string(this->ArgV[this->Index] ) ) );
	}
      item->Evaluate( this->ArgC, this->ArgV, this->Index );
      } 
    else
      {
      const char* optChar = this->ArgV[this->Index]+1;
      while ( *optChar ) 
	{
	// short option
	for ( KeyActionListType::iterator it = this->m_KeyActionListComplete.begin(); it != this->m_KeyActionListComplete.end(); ++it )
	  {
	  if ( (*it)->m_Key == *optChar )
	    {
	    item = (*it)->m_Action;
	    }
	  }
      
	if ( !item ) 
	  {
	  const char opt[2] = { *optChar, 0 };
	  throw( Exception( std::string("Unknown option: -") + std::string(opt) ) );
	  }
	item->Evaluate( this->ArgC, this->ArgV, this->Index );
	
	++optChar; // next short option in this block, POSIX style.
	}
      }
    
    ++this->Index;
    } // while

  for ( NonOptionParameterListType::iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it, ++this->Index )
    {
    if ( this->Index >= this->ArgC )
      {
      throw( Exception( "Insufficient number of command line arguments", this->Index ) );
      }
    (*it)->Evaluate( this->ArgC, this->ArgV, this->Index );
    }
  
  return true;
}

void
CommandLine::PrintHelp
() const
{
  ProgramPropertiesMapType::const_iterator it = this->m_ProgramInfo.find(PRG_TITLE);
  if ( it != this->m_ProgramInfo.end() )
    {
    StdErr << "TILE:\n\n";
    StdErr.FormatText( it->second, 5 ) << "\n";
    }

  it = this->m_ProgramInfo.find(PRG_DESCR);
  if ( it != this->m_ProgramInfo.end() )
    {
    StdErr << "DESCRIPTION:\n\n";
    StdErr.FormatText( it->second, 5 ) << "\n";
    }

  it = this->m_ProgramInfo.find(PRG_SYNTX);
  if ( it != this->m_ProgramInfo.end() )
    {
    StdErr << "SYNTAX:\n\n";
    StdErr.FormatText( it->second, 5 ) << "\n";
    }

  StdErr << "LIST OF SUPPORTED OPTIONS:\n\n";

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

Console& operator<<( Console& console, CommandLine::Exception e )
{
  console << e.Message << " [argument #" << e.Index << "]\n";
  return console;
}

} // namespace cmtk
