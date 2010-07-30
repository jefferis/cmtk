/*
//
//  Copyright 2010 SRI International
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

#include "cmtkSimpleLevelsetDevice_kernels.h"

#include "GPU/cmtkCUDA.h"

void
cmtk::SimpleLevelsetDeviceUpdateInsideOutside( float* levelset, float* volume, const int nPixels, float* insideSum, float* outsideSum, int* nInside )
{
}

void
cmtk::SimpleLevelsetDeviceUpdateLevelset( float* levelset, float* volume, const int nPixels, const float mInside, const float mOutside, const float timeDelta, const float levelsetThreshold )
{
}
