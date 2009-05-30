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

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class T,ScalarDataType DT,cmtk::Interpolators::InterpolationEnum I>
T inline VoxelMatchingMetric<T,DT,I>::GetSampleY
( const size_t baseIndex, const Types::Coordinate* frac ) 
  const
{
  const Types::Coordinate offsX = 1.0 - frac[0];
  const Types::Coordinate offsY = 1.0 - frac[1];
  const Types::Coordinate offsZ = 1.0 - frac[2];
  
  assert( (baseIndex+this->DataY.nextIJK) < this->DataY.NumberOfSamples );
  const T *node = this->DataY.Data + baseIndex;
  
  return static_cast<T>
    ( offsZ*(offsY*(offsX*node[0] + frac[0]*node[1]) +
	     frac[1]*(offsX*node[this->DataY.nextJ] + 
		      frac[0]*node[this->DataY.nextIJ]) )+
      frac[2]*(offsY*(offsX*node[this->DataY.nextK] + 
		      frac[0]*node[this->DataY.nextIK]) +
	       frac[1]*(offsX*node[this->DataY.nextJK] + 
			frac[0]*node[this->DataY.nextIJK]) ) );
}

} // namespace cmtk
