/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include "cmtkLandmarkPair.h"

std::ostream& 
operator<<( std::ostream& stream, const cmtk::LandmarkPair& pair )
{
  stream << pair.m_Location << "\t" << pair.m_TargetLocation << "\t" << pair.m_Name << std::endl;
  return stream;
}

std::istream& 
operator>>( std::istream& stream, cmtk::LandmarkPair& pair )
{
  stream >> pair.m_Location >> pair.m_TargetLocation >> pair.m_Name;
  return stream;
}

