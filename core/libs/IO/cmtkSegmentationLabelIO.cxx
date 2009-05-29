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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkSegmentationLabelIO.h>

#include <string>
#include <sstream>

namespace
cmtk
{

/** \addtogroup IO */
//@{

std::istream& operator>>
( std::istream& stream, SegmentationLabelMap& lblMap )
{
  std::string line;
  
  while ( ! stream.eof() ) 
    {
    std::getline( stream, line );
    if ( line.length() && (line[0] != '#') ) 
      { // skip blank and comments
      int id;
      std::string name, rs, gs, bs, as;
      
      std::istringstream inStr( line );
      inStr >> id >> name >> rs >> gs >> bs >> as;
      
      lblMap[id].SetName( name.c_str() );
      lblMap[id].SetRGB( atoi( rs.c_str() ), atoi( gs.c_str() ), atoi( bs.c_str() ) );
      }
    }
  
  return stream;
}

} // namespace cmtk
