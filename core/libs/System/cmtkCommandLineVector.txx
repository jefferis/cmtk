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
#include <sstream>

template<class T>
void
cmtk::CommandLine::Vector<T>
::Evaluate( const size_t argc, const char* argv[], size_t& index )
{
  if ( !this->m_HasBeenUsed )
    {
    this->m_pVector->resize(0);
    this->m_HasBeenUsed = true;
    }

  if ( index+1 < argc ) 
    {
    ++index;
    // first, replace all commas with spaces, so we can simply use a stringstream for parsing the vector elements
    std::string str( argv[index] );
    for ( int i = 0; i < str.length(); ++i )
      {
      if ( str[i] == ',' )
	str[i] = ' ';
      }

    // new read values from the whitespaced argument
    std::istringstream strm( str );
    while ( strm.good() && ! strm.eof() )
      {
      T nextValue;
      strm >> nextValue;
      this->m_pVector->push_back( nextValue );
      }

    } 
  else
    {
    throw( Exception( "Vector command line option needs an argument.", index ) );
    }
}

template<class T>
mxml_node_t* 
cmtk::CommandLine::Vector<T>
::MakeXML(  mxml_node_t *const parent ) const 
{
  if ( ! (this->m_Properties & PROPS_NOXML) )
    {
    const std::string typeName = std::string ( CommandLineTypeTraits<T>::GetName() ) + std::string( "-vector" );    
    mxml_node_t *node = mxmlNewElement( parent, typeName.c_str() );
    
    // write any attributes the user might have set
    for ( std::map<const std::string,std::string>::const_iterator attrIt = this->m_Attributes.begin(); attrIt != this->m_Attributes.end(); ++attrIt )
      {
      mxmlElementSetAttr( node, attrIt->first.c_str(), attrIt->second.c_str() );
      }
    
    mxmlElementSetAttr( node, "multiple", "true" );
    
    return node;
    }
  return NULL;
}
