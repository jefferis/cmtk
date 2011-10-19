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
CommandLine::PrintHelp
( const bool advanced ) const
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
      if ( (((*grp)->GetProperties() & Self::PROPS_ADVANCED)==0) || advanced )
	{
	StdOut << (*grp)->m_Description << "\n\n";
	
	const size_t indent = 2;
	
	const KeyActionListType& kal = (*grp)->m_KeyActionList;
	for ( KeyActionListType::const_iterator it = kal.begin(); it != kal.end(); ++it )
	  {
	  (*it)->PrintHelp( indent, advanced );
	  }
	}
      }
    }
  
  StdOut << "\n";
}

} // namespace cmtk
