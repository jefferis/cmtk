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
CommandLine::PrintMan
( const char* argv0 ) const
{
  // determine program file name
  const char* cmd = strrchr( argv0, '/' );
  if ( cmd )
    {
    ++cmd;
    }
  else
    {
    cmd = argv0;
    }

  ProgramPropertiesMapType::const_iterator ppit = this->m_ProgramInfo.find(PRG_VERSN);
  StdOut << ".TH " << cmd << " \"1\" \"" << __DATE__ << "\" \"CMTK " << ppit->second << "\" \"The Computational Morphometry Toolkit\"\n";

  ppit = this->m_ProgramInfo.find(PRG_TITLE);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << ".SH NAME\n" << cmd << " \\- " << ppit->second << "\n";
    }

  ppit = this->m_ProgramInfo.find(PRG_SYNTX);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << ".SH SYNOPSIS\n";
    StdOut << ppit->second << "\n";
    }
  else
    {
    if ( this->m_NonOptionParameterList.size() || this->m_NonOptionParameterVectorList.size() )
      {
      StdOut << ".SH SYNOPSIS\n\\fB" << cmd << "\\fR ";
      
      for ( NonOptionParameterListType::const_iterator it = this->m_NonOptionParameterList.begin(); it != this->m_NonOptionParameterList.end(); ++it )
	{
	StdOut << (*it)->m_Name << " ";
	}
      for ( NonOptionParameterVectorListType::const_iterator it = this->m_NonOptionParameterVectorList.begin(); it != this->m_NonOptionParameterVectorList.end(); ++it )
	{
	StdOut << (*it)->m_Name << " ";
	}
      StdOut << "\n";
      }
    }

  ppit = this->m_ProgramInfo.find(PRG_DESCR);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << ".SH DESCRIPTION\n";
    StdOut << ppit->second << "\n";
    }
  
  StdOut << ".SH OPTIONS\n";

  for ( KeyActionGroupListType::const_iterator grp = this->m_KeyActionGroupList.begin(); grp != this->m_KeyActionGroupList.end(); ++grp )
    {
    if ( ! (*grp)->m_KeyActionList.empty() )
      {
      StdOut << ".SS " << (*grp)->m_Description << "\n";

      const KeyActionListType& kal = (*grp)->m_KeyActionList;
      for ( KeyActionListType::const_iterator it = kal.begin(); it != kal.end(); ++it )
	{
	(*it)->PrintManWithPrefix();
	}
      }
    }

  ppit = this->m_ProgramInfo.find(PRG_CNTRB);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << ".SH AUTHORS\n" << ppit->second << "\n";
    }

  ppit = this->m_ProgramInfo.find(PRG_LCNSE);
  if ( ppit != this->m_ProgramInfo.end() )
    {
    StdOut << ".SH LICENSE\n" << ppit->second << "\n";
    }

  StdOut << ".SH BUGS\nReport bugs at http://nitrc.org/projects/cmtk/\n";
  StdOut << ".SH ACKNOWLEDGMENTS\n"
	 << "From April 2009 through September 2011, CMTK Development and Maintenance was supported by the National Institute of Biomedical Imaging and Bioengineering under Grant No.R01 EB008381 (PI: Torsten Rohlfing).\n";
}

} // namespace cmtk
