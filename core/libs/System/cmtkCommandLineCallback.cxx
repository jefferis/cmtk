/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2009, 2013 SRI International
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

#include <System/cmtkCoverity.h>

void 
cmtk::CommandLine::Callback::Evaluate
( const size_t argc, const char* argv[], size_t& index )
{
  // callback with argument?
  if ( this->m_FuncArg ) 
    {
    if ( index+1 < argc ) 
      {
      try
	{
	this->m_FuncArg( argv[index+1] );
	}
      catch ( const char* error )
	{
	throw( Exception( error, index ) );
	} 
      ++index;
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
	try
	  {
	  this->m_FuncIntArg( ConvertStrToLong( argv[index+1] ) );
	  }
	catch ( const char* error )
	  {
	  throw( Exception( error, index ) );
	  } 

	++index;
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
	try
	  {
	  this->m_FuncDblArg( ConvertStrToDouble( argv[index+1] ) );
	  }
	catch ( const char* error )
	  {
	  throw( Exception( error, index ) );
	  } 

	++index;
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
	try
	  {
	  this->m_FuncMultiArg( argv+index+1, argsUsed );
	  }
	catch ( const char* error )
	  {
	  throw( Exception( error, index ) );
	  } 
	index += argsUsed;
	} 
      else
	{
	throw( Exception( "Option needs an argument", index ) );
	}
      } 
    else
      {
      // no argument to callback
      try
	{
	this->m_Func();
	}
      catch ( const char* error )
	{
	throw( Exception( error, index ) );
	}
      }
}

mxml_node_t* 
cmtk::CommandLine::Callback
::MakeXML(  mxml_node_t *const parent ) const 
{
  mxml_node_t* node = NULL;
  if ( this->m_Func )
    {
    node = mxmlNewElement( parent, "boolean" );
    mxml_node_t *dflt = mxmlNewElement( node, "default" );
    Coverity::FakeFree( mxmlNewText( dflt, 0, "false" ) );
    }
  else if ( this->m_FuncArg )
    {
    node = mxmlNewElement( parent, "string" );
    }
  else if ( this->m_FuncIntArg )
    {
    node = mxmlNewElement( parent, "integer" );
    }
  else if ( this->m_FuncDblArg )
    {
    node = mxmlNewElement( parent, "double" );
    }
  else if ( this->m_FuncMultiArg )
    {
    node = mxmlNewElement( parent, "string-vector" );
    }

  mxmlElementSetAttr( node, "multiple", "true" );
  return node;
}

std::string
cmtk::CommandLine::Callback
::GetParamTypeString() const
{
  if ( this->m_FuncArg )
    {
    return Item::Helper<std::string>::GetParamTypeString( this );
    }
  else if ( this->m_FuncIntArg )
    {
    return Item::Helper<int>::GetParamTypeString( this );
    }
  else if ( this->m_FuncDblArg )
    {
    return Item::Helper<double>::GetParamTypeString( this );
    }
  else if ( this->m_FuncMultiArg )
    {
    return "<string-vector>";
    }
  return "";
}
