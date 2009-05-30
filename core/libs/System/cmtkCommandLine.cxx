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
}

void
CommandLine
::BeginGroup( const char* name, const char* description ) 
{ 
  this->m_KeyActionGroupList.push_back( KeyActionGroupType::SmartPtr( new KeyActionGroupType( name, description ) ) );
  this->m_KeyActionList = &(this->m_KeyActionGroupList.back()->m_KeyActionList);
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
CommandLine::Parse()
{
  Index = 1;
  while ( (Index < ArgC) && (ArgV[Index][0] == '-') ) 
    {
    // Break at first non-switch argument.
    if ( ArgV[Index][0] != '-' ) return true;
    
    // Check for "--xml" special option, which produces self description according to Slicer execution model.
    if ( !strcmp( ArgV[Index], "--xml" ) ) 
      {
      this->WriteXML();
      exit( 2 );
      }

    // Check for "--help" special option, which produces textual description of all command line options
    if ( !strcmp( ArgV[Index], "--help" ) ) 
      {
      this->PrintHelp();
      exit( 2 );
      }

    // Like POSIX, break at "--" terminator.
    if ( !strcmp( ArgV[Index], "--" ) ) 
      {
      ++Index;
      break;
      }
    
    SmartPointer<Item> item(NULL);
    if ( ArgV[Index][1] == '-' ) 
      {
      // long option
      for ( KeyActionListType::iterator it = this->m_KeyActionListComplete.begin(); it != this->m_KeyActionListComplete.end(); ++it )
	{
	if ( (*it)->m_KeyString == std::string( ArgV[Index]+2 ) )
	  {
	  item = (*it)->m_Action;
	  }
	}
      
      // not found?
      if ( !item ) 
	{
	StdErr.printf( "Unknown option: %s\n", ArgV[Index] );
	return false;
	}
      item->Evaluate( ArgC, ArgV, Index );
      } 
    else
      {
      const char* optChar = ArgV[Index]+1;
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
	  StdErr.printf( "Unknown option: -%c\n", *optChar );
	  return false;
	  }
	item->Evaluate( ArgC, ArgV, Index );
	
	++optChar; // next short option in this block, POSIX style.
	}
      }
    
    ++Index;
    } // while
  
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
