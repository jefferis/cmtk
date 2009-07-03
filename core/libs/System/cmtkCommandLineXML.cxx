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
#include <cmtkconfig.h>

#include <cmtkCommandLine.h>

#include <mxml.h>

#include <stdio.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

mxml_node_t*
CommandLine::AddProgramInfoXML( mxml_node_t *const parent, const ProgramProperties key, const char* name ) const
{
  ProgramPropertiesMapType::const_iterator it = this->m_ProgramInfo.find( key );
  if ( it != this->m_ProgramInfo.end() )
    {
    mxml_node_t *node = mxmlNewElement( parent, name );
    mxmlNewText( node, 0, it->second.c_str() );
    return node;
    }
  return NULL;
}

const char *
cmtkWhitespaceWriteMiniXML( mxml_node_t*, int where)
{
  switch ( where )
    {
    case MXML_WS_BEFORE_OPEN:
      return NULL;
    case MXML_WS_AFTER_OPEN:
      return "\n";
    case MXML_WS_BEFORE_CLOSE:
      return "\n";
    case MXML_WS_AFTER_CLOSE:
      return "\n";
    }
  return NULL;
}

void
CommandLine::WriteXML
() const
{
  // check for globally disabled XML support (some tools should not be used as Slicer plugins, for example)
  if ( ! (this->m_Properties & PROPS_NOXML) )
    {
    mxml_node_t *xml = mxmlNewXML("1.0");
    
    mxml_node_t *x_exec = mxmlNewElement(xml, "executable");
    
    this->AddProgramInfoXML( x_exec, PRG_CATEG, "category" );
    this->AddProgramInfoXML( x_exec, PRG_TITLE, "title" );
    this->AddProgramInfoXML( x_exec, PRG_DESCR, "description" );
    this->AddProgramInfoXML( x_exec, PRG_LCNSE, "license" );
    this->AddProgramInfoXML( x_exec, PRG_CNTRB, "contributor" );
    this->AddProgramInfoXML( x_exec, PRG_ACKNL, "acknowledgements" );
    this->AddProgramInfoXML( x_exec, PRG_DOCUM, "documentation-url" );
    this->AddProgramInfoXML( x_exec, PRG_VERSN, "version" );
    
    for ( KeyActionGroupListType::const_iterator grp = this->m_KeyActionGroupList.begin(); grp != this->m_KeyActionGroupList.end(); ++grp )
      {
      if ( ! ((*grp)->GetProperties() & PROPS_NOXML) )
	{
	mxml_node_t *parameterGroup = mxmlNewElement( x_exec, "parameters" );
	
	if ( (*grp)->GetProperties() & PROPS_ADVANCED )
	  mxmlElementSetAttr( parameterGroup, "advanced", "true" );
	
	const std::string& name = (*grp)->m_Name;
	if ( name == "MAIN" )
	  {
	  mxmlNewText( mxmlNewElement( parameterGroup, "label" ), 0, "General" );
	  mxmlNewText( mxmlNewElement( parameterGroup, "description" ), 0, "General Parameters" );
	  
	  int index = 0;
	  for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	    {
	    (*it)->MakeXML( parameterGroup, index++ );
	    }
	  }
	else
	  {
	  mxmlNewText( mxmlNewElement( parameterGroup, "label" ), 0, name.c_str() );
	  mxmlNewText( mxmlNewElement( parameterGroup, "description" ), 0, (*grp)->m_Description.c_str() );
	  }
	
	const KeyActionListType& kal = (*grp)->m_KeyActionList;
	for ( KeyActionListType::const_iterator it = kal.begin(); it != kal.end(); ++it )
	  {
	  (*it)->MakeXML( parameterGroup );
	  }
	}
      }
    
    mxmlSaveFile( xml, stdout, cmtkWhitespaceWriteMiniXML );
    fputs( "\n", stdout ); // Slicer's XML parser needs an extra \n after the last line

    mxmlDelete( xml );
    }
}

} // namespace cmtk
