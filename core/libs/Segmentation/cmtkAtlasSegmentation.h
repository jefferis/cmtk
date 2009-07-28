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

#include <cmtkUniformVolume.h>

namespace
cmtk
{

class AtlasSegmentation
{
public:
  /// Constructor: compute registrations.
  AtlasSegmentation( UniformVolume::SmartPtr& targetImage, UniformVolume::SmartPtr& atlasImg, UniformVolume::SmartPtr& atlasLbl, const bool verbose );

  /// Get reformatted label map.
  cmtk::UniformVolume::SmartPtr& GetLabelMap()
  {
    return this->m_LabelMap;
  }

private:
  /// Output label map.
  UniformVolume::SmartPtr m_LabelMap;
};

} // namespace cmtk
