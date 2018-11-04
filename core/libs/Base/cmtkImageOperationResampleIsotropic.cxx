/*
//
//  Copyright 2009-2011, 2013 SRI International
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

#include "cmtkImageOperationResampleIsotropic.h"

cmtk::UniformVolume::SmartPtr cmtk::ImageOperationResampleIsotropic ::Apply(
    cmtk::UniformVolume::SmartPtr &volume) {
  if (this->m_Exact)
    return cmtk::UniformVolume::SmartPtr(
        volume->GetResampledExact(this->m_Resolution));
  else
    return cmtk::UniformVolume::SmartPtr(
        volume->GetResampled(this->m_Resolution, true /*allowUpsampling*/));
}

void cmtk::ImageOperationResampleIsotropic ::New(const double resolution) {
  ImageOperation::m_ImageOperationList.push_back(
      SmartPtr(new ImageOperationResampleIsotropic(resolution)));
}

void cmtk::ImageOperationResampleIsotropic ::NewExact(const double resolution) {
  ImageOperation::m_ImageOperationList.push_back(
      SmartPtr(new ImageOperationResampleIsotropic(resolution, true)));
}
