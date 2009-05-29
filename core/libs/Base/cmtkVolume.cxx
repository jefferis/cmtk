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

#include <cmtkVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

bool
Volume::GetTrilinear
( ProbeInfo& probeInfo, const int X, const int Y, const int Z,
  const Vector3D& Location, const Types::Coordinate* from, 
  const Types::Coordinate* to ) const
{
  const TypedArray* data = this->GetData();

  int offset = X+Dims[0]*(Y+Dims[1]*Z);

  bool data_present = data->Get( probeInfo.Values[0], offset );
  
  if ( X<Dims[0]-1 ) 
    {
    data_present &= data->Get( probeInfo.Values[1], offset+nextI );
    
    if ( Y<Dims[1]-1 ) 
      {
      data_present &= data->Get( probeInfo.Values[3], offset+nextIJ );
      
      if ( Z<Dims[2]-1 )
	data_present &= data->Get( probeInfo.Values[7], offset+nextIJK );
      }
    if ( Z<Dims[2]-1 )
      data_present &= data->Get( probeInfo.Values[5], offset+nextIK );
    }
  
  if ( Y<Dims[1]-1 ) 
    {
    data_present &= data->Get( probeInfo.Values[2], offset+nextJ );
    
    if ( Z<Dims[2]-1 )
      data_present &= data->Get( probeInfo.Values[6], offset+nextJK );
    }
  
  if ( Z<Dims[2]-1 )
    data_present &= data->Get( probeInfo.Values[4], offset+nextK );
  
  if (data_present)
    {
    for ( int i=0; i<3; ++i ) 
      {
      probeInfo.Deltas[i] = 1.0/(to[i]-from[i]);
      
      probeInfo.Offsets[i] = 1- (probeInfo.Offsets[3+i] = probeInfo.Deltas[i]*(Location[i]-from[i]) );
      }
    
    probeInfo.Location = Location;
    
    return true;
    }
  
  return false;
}

Vector3D
Volume::GetCenter () const 
{
  return this->m_Origin + 0.5 * Vector3D(Size);
}

} // namespace cmtk
