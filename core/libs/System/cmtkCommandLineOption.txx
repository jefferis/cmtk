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

#include <typeinfo>

template<class T>
void
cmtk::CommandLine::Option<T>
::Evaluate( const size_t argc, const char* argv[], size_t& index )
{
  if ( this->Flag ) 
    *(this->Flag) = true;
  if ( index+1 < argc ) 
    {
    *(this->Var) = this->Convert<T>( argv[index+1] );
	++index;
    } 
  else
    {
    throw( Exception( "Option needs an argument.", index ) );
    }
}

template<class T>
mxml_node_t* 
cmtk::CommandLine::Option<T>
::MakeXML(  mxml_node_t *const parent ) const 
{
  if ( Flag ) // if there is a flag monitoring this option, then we need to create an additional boolean toggle
    {
    }
  
  const char* typeName = CommandLineTypeTraits<T>::GetName();

  mxml_node_t *node = NULL;
  if ( std::string( typeName ) == "string" )
    {
    if ( this->m_Properties & PROPS_IMAGE )
      {
      node = mxmlNewElement( parent, "image" );
      
      if ( this->m_Properties & PROPS_LABELS )
	mxmlElementSetAttr( node, "type", "label" );
      else
	mxmlElementSetAttr( node, "type", "scalar" );
      }
    else if ( this->m_Properties & PROPS_XFORM )
      {
      node = mxmlNewElement( parent, "transform" );
      mxmlElementSetAttr( node, "fileExtensions", ".txt" );
      }
    else if ( this->m_Properties & PROPS_FILENAME )
      node = mxmlNewElement( parent, "file" );
    else if ( this->m_Properties & PROPS_DIRNAME )
      node = mxmlNewElement( parent, "directory" );
    else 
      node = mxmlNewElement( parent, "string" );

    if ( this->m_Properties & PROPS_OUTPUT )
      mxmlNewText( mxmlNewElement( node, "channel" ), 0, "output" );
    else
      mxmlNewText( mxmlNewElement( node, "channel" ), 0, "input" );
    }
  else
    node = mxmlNewElement( parent, typeName );

  // write any attributes the user might have set
  for ( std::map<const std::string,std::string>::const_iterator attrIt = this->m_Attributes.begin(); attrIt != this->m_Attributes.end(); ++attrIt )
    {
    mxmlElementSetAttr( node, attrIt->first.c_str(), attrIt->second.c_str() );
    }
  
  if ( !Flag ) // if there is no flag monitoring this option, then there must be a valid default value
    {
    mxml_node_t *dflt = mxmlNewElement( node, "default" );
    mxmlNewText( dflt, 0, CommandLineTypeTraits<T>::ValueToStringMinimal( this->Var ).c_str() );
    }
  return node;
}

template<class T>
std::ostringstream& 
cmtk::CommandLine::Option<T>
::PrintHelp( std::ostringstream& fmt ) const
{
  if ( this->Flag && !(*this->Flag) )
    fmt << "\n[Default: disabled]";
  else
    fmt << "\n[Default value: " << CommandLineTypeTraits<T>::ValueToString( this->Var ) << "]";
  return fmt;
}

template<class T>
void
cmtk::CommandLine::Option<T>
::PrintWiki() const
{
  if ( this->Flag && !(*this->Flag) )
    StdOut << " '''[Default: disabled]'''";
  else
    StdOut << " '''[Default value: " << CommandLineTypeTraits<T>::ValueToString( this->Var ) << "]'''";
}
