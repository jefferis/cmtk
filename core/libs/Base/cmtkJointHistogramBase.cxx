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

#include "cmtkJointHistogramBase.h"

size_t
cmtk::JointHistogramBase::CalcNumBins
( const size_t numberOfSamples, const Types::DataItemRange& valueRange ) 
{
  const size_t side = static_cast<size_t>( sqrt( static_cast<float>( numberOfSamples )) );
  const size_t dataRange = static_cast<size_t>( valueRange.Width() + 1 );
  const int upperLimit = std::min<size_t>( std::min<size_t>( dataRange, 128), side );
  return std::max<size_t>( 8, upperLimit );
}

size_t 
cmtk::JointHistogramBase::CalcNumBins ( const UniformVolume* volume ) 
{
  return Self::CalcNumBins( volume->CropRegion().Size(), volume->GetData()->GetRange() ) ;
}
