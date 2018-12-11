/*
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

#include "cmtkVoxelMatchingMetric.h"

#include <Base/cmtkJointHistogram.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class T,ScalarDataType DT,Interpolators::InterpolationEnum I>
VoxelMatchingMetric<T,DT,I>::VoxelMatchingMetric
( const UniformVolume* refVolume, const UniformVolume* fltVolume, const bool initData )
{
  this->DataX.PrecomputeIncrements( refVolume );
  this->DataY.PrecomputeIncrements( fltVolume );

  if ( initData ) 
    {
    this->DataX.Init( refVolume );
    this->DataY.Init( fltVolume );
    }
}

// instantiate necessary templates.
template class VoxelMatchingMetric<short,TYPE_SHORT,Interpolators::LINEAR>;
template class VoxelMatchingMetric<byte,TYPE_BYTE,Interpolators::LINEAR>;

template class VoxelMatchingMetric<short,TYPE_SHORT,Interpolators::NEAREST_NEIGHBOR>;
template class VoxelMatchingMetric<byte,TYPE_BYTE,Interpolators::NEAREST_NEIGHBOR>;

} // namespace cmtk
