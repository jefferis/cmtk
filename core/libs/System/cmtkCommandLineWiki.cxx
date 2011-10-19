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
#include <cmtkconfig.h>

#include "cmtkCommandLine.h"

#include <System/cmtkConsole.h>

namespace
cmtk
{

/** \addtogroup System */
//@{

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

} // namespace cmtk
