/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2010 SRI International
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

#include "cmtkUniformVolumeInterpolatorBase.h"

#include <limits>

void cmtk::UniformVolumeInterpolatorBase ::SetVolume(
    const UniformVolume &volume) {
  const TypedArray &data = *(volume.GetData());
  const size_t nPixels = data.GetDataSize();
  this->m_VolumeDataArray.resize(nPixels);
  for (Types::GridIndexType n = 0; n < nPixels; ++n) {
    if (!data.Get(this->m_VolumeDataArray[n], n))
      this->m_VolumeDataArray[n] =
          std::numeric_limits<Types::DataItem>::infinity();
  }

  this->m_VolumeDims = volume.GetDims();
  this->m_VolumeDeltas = volume.Deltas();
  this->m_VolumeOffset = volume.m_Offset;
  this->m_NextJ = this->m_VolumeDims[0];
  this->m_NextK = this->m_NextJ * this->m_VolumeDims[1];
}
