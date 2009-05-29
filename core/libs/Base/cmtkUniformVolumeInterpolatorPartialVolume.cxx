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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

#include <cmtkUniformVolumeInterpolatorPartialVolume.h>

namespace cmtk
{

/** \addtogroup Base */
//@{

bool
UniformVolumeInterpolatorPartialVolume
::GetDataAt(const Vector3D& v, Types::DataItem& value) const
{
  value=0;

  const int *imageDims = this->m_Volume->GetDims();

  const Types::Coordinate *Delta = this->m_Volume->GetDelta();

  Types::Coordinate lScaled[3];
  int imageGridPoint[3];
  for ( int n = 0; n < 3; ++n )
    {
    lScaled[n] = (v[n]-this->m_Volume->m_Origin[n]) / Delta[n];
    imageGridPoint[n] = (int) floor( lScaled[n] );
    if ( ( imageGridPoint[n] < 0 ) || ( imageGridPoint[n] >= imageDims[n]-1 ) )
      return false;
    }
  
  const size_t offset = imageGridPoint[0] + imageDims[0] * ( imageGridPoint[1] + imageDims[1] * imageGridPoint[2]);
  const TypedArray* gridData = this->m_Volume->GetData();
  
  Types::DataItem corners[8];
  bool data_present = gridData->Get( corners[0], offset );
  data_present &= gridData->Get( corners[1], offset+this->m_Volume->GetNextI() );
  data_present &= gridData->Get( corners[2], offset+this->m_Volume->GetNextJ() );
  data_present &= gridData->Get( corners[3], offset+this->m_Volume->GetNextIJ() );
  data_present &= gridData->Get( corners[4], offset+this->m_Volume->GetNextK() );
  data_present &= gridData->Get( corners[5], offset+this->m_Volume->GetNextIK() );
  data_present &= gridData->Get( corners[6], offset+this->m_Volume->GetNextJK() );
  data_present &= gridData->Get( corners[7], offset+this->m_Volume->GetNextIJK() );
  
  if (data_present) 
    {
    const Types::Coordinate revX = lScaled[0]-imageGridPoint[0];
    const Types::Coordinate revY = lScaled[1]-imageGridPoint[1];
    const Types::Coordinate revZ = lScaled[2]-imageGridPoint[2];
    const Types::Coordinate offsX = 1-revX;
    const Types::Coordinate offsY = 1-revY;
    const Types::Coordinate offsZ = 1-revZ;

    const Types::Coordinate weights[8] = 
      { offsX * offsY * offsZ, revX * offsY * offsZ, offsX * revY * offsZ, revX * revY * offsZ,
	offsX * offsY * revZ, revX * offsY * revZ, offsX * revY * revZ, revX * revY * revZ 
      };

    bool done[8];
    memset( done, 0, sizeof( done ) );
    
    Types::Coordinate maxWeight = 0;
    for ( unsigned int j = 0; j < 8; ++j ) 
      {
      if ( done[j] ) continue;
      Types::Coordinate weight = weights[j];
      for ( unsigned int i = j+1; i < 8; ++i ) 
	{
	if ( done[i] ) continue;
	if ( corners[i] == corners[j] ) 
	  {
	  weight += weights[i];
	  done[i] = true;
	  }
	}
      if ( weight > maxWeight ) 
	{
	value = corners[j];
	maxWeight = weight;
	}
      }
    
    return true;
    }
  
  return false;
}

} // namespace cmtk

