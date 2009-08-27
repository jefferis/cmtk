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
mxml_node_t* 
cmtk::CommandLine::Item::Helper<T>
::MakeXML( const Item* item, mxml_node_t *const parent )
{
  if ( ! (item->m_Properties & PROPS_NOXML) )
    {
    const char* typeName = CommandLineTypeTraits<T>::GetName();
    
    mxml_node_t *node = NULL;
    if ( std::string( typeName ) == "string" )
      {
      if ( item->m_Properties & PROPS_IMAGE )
	{
	node = mxmlNewElement( parent, "image" );
	
	if ( item->m_Properties & PROPS_LABELS )
	  mxmlElementSetAttr( node, "type", "label" );
	else
	  mxmlElementSetAttr( node, "type", "scalar" );
	}
      else if ( item->m_Properties & PROPS_XFORM )
	{
	node = mxmlNewElement( parent, "transform" );
	mxmlElementSetAttr( node, "fileExtensions", ".txt" );
	}
      else if ( item->m_Properties & PROPS_FILENAME )
	node = mxmlNewElement( parent, "file" );
      else if ( item->m_Properties & PROPS_DIRNAME )
	node = mxmlNewElement( parent, "directory" );
      else 
	node = mxmlNewElement( parent, "string" );
      
      if ( item->m_Properties & PROPS_OUTPUT )
	mxmlNewText( mxmlNewElement( node, "channel" ), 0, "output" );
      else
	mxmlNewText( mxmlNewElement( node, "channel" ), 0, "input" );
      }
    else
      node = mxmlNewElement( parent, typeName );
    
    // write any attributes the user might have set
    for ( std::map<const std::string,std::string>::const_iterator attrIt = item->m_Attributes.begin(); attrIt != item->m_Attributes.end(); ++attrIt )
      {
      mxmlElementSetAttr( node, attrIt->first.c_str(), attrIt->second.c_str() );
      }
    
    return node;
    }
  return NULL;
}

template<class T>
std::string
cmtk::CommandLine::Item::Helper<T>
::GetParamTypeString( const Item* item )
{
  const std::string& typeName = CommandLineTypeTraits<T>::GetName();
    
  if ( typeName == "string" )
    {
    if ( item->m_Properties & PROPS_IMAGE )
      {
      if ( item->m_Properties & PROPS_LABELS )
	return "<labelmap-path>";
      else
	return "<image-path>";
      }
    else if ( item->m_Properties & PROPS_XFORM )
      {
      return "<transformation-path>";
      }
    else if ( item->m_Properties & PROPS_FILENAME )
      return "<path>";
    else if ( item->m_Properties & PROPS_DIRNAME )
      return "<directory>";
    else 
      return "<string>";
    }

  return std::string("<")+typeName+std::string(">");
}
