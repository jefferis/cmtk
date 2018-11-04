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

#include "cmtkLandmarkList.h"

namespace cmtk {

/** \addtogroup Base */
//@{

LandmarkList::Iterator LandmarkList::FindByName(const std::string &name) {
  Self::Iterator it = this->begin();
  while (it != this->end()) {
    if (it->m_Name == name) {
      return it;
    }
    ++it;
  }

  return it;
}

LandmarkList::ConstIterator LandmarkList::FindByName(
    const std::string &name) const {
  Self::ConstIterator it = this->begin();
  while (it != this->end()) {
    if (it->m_Name == name) {
      return it;
    }
    ++it;
  }

  return it;
}

}  // namespace cmtk
