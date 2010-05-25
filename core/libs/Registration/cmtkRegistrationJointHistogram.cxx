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

#include <cmtkRegistrationJointHistogram.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<cmtk::Interpolators::InterpolationEnum I>
RegistrationJointHistogram<I>::RegistrationJointHistogram 
( const UniformVolume* refVolume, const UniformVolume* fltVolume,
  const unsigned int numBinsX, const unsigned int numBinsY,
  const Types::DataItemRange& boundsX,
  const Types::DataItemRange& boundsY ) :
#ifdef CMTK_PVI_HISTOGRAMS
  JointHistogram<float>(),
#else
  JointHistogram<int>(),
#endif
  VoxelMatchingMetric<byte,TYPE_BYTE,I>( refVolume, fltVolume, false /* initData */ )
{
  this->SetNumBins( this->DataX.Init( refVolume, numBinsX, boundsX ), this->DataY.Init( fltVolume, numBinsY, boundsY ) );
}

template<cmtk::Interpolators::InterpolationEnum I>
unsigned int
RegistrationJointHistogram<I>::SetDataX 
( const UniformVolume* volume, const unsigned int numBins, 
  const Types::DataItemRange& bounds )
{
  this->VoxelMatchingMetric<byte,TYPE_BYTE,I>::SetDataX( volume );
  this->SetNumBinsX( this->DataX.Init( volume, numBins, bounds ) );
  return NumBinsX;
}

template<cmtk::Interpolators::InterpolationEnum I>
unsigned int
RegistrationJointHistogram<I>::SetDataY 
( const UniformVolume* volume, const unsigned int numBins, 
  const Types::DataItemRange& bounds )
{
  this->VoxelMatchingMetric<byte,TYPE_BYTE,I>::SetDataY( volume );
  this->SetNumBinsY( this->DataY.Init( volume, numBins, bounds ) );
  return NumBinsY;
}

template<cmtk::Interpolators::InterpolationEnum I>
void
RegistrationJointHistogram<I>::SetDataXY
( const UniformVolume* volumeX, const unsigned int numBinsX, 
  const UniformVolume* volumeY, const unsigned int numBinsY,
  const Types::DataItemRange& boundsX, const Types::DataItemRange& boundsY )
{
  this->VoxelMatchingMetric<byte,TYPE_BYTE,I>::SetDataXY( volumeX, volumeY );
  this->SetNumBins( this->DataX.Init( volumeX, numBinsX, boundsX ), this->DataY.Init( volumeY, numBinsY, boundsY ) );
}

// instantiate required templates
template class RegistrationJointHistogram<Interpolators::LINEAR>;
template class RegistrationJointHistogram<Interpolators::NEAREST_NEIGHBOR>;

} // namespace cmtk
