/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2012 SRI International
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

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class TData> 
inline bool
UniformVolume::ProbeData 
( TData& result, const TData* dataPtr, const Self::CoordinateVectorType& location ) const
{
  result=0;
  
  Self::CoordinateVectorType l( location );
  l -= this->m_Offset;

  if ( (l[0] < 0) || (l[1] < 0) || (l[2] < 0) )
    return false;
  
  const Types::GridIndexType idxX=(Types::GridIndexType) (l[0]/this->m_Delta[0]);
  if ( idxX>=this->m_Dims[0]-1 )
    return false;
  
  const Types::GridIndexType idxY=(Types::GridIndexType) (l[1]/this->m_Delta[1]);
  if ( idxY>=this->m_Dims[1]-1 )
    return false;

  const Types::GridIndexType idxZ=(Types::GridIndexType) (l[2]/this->m_Delta[2]);
  if ( idxZ>=this->m_Dims[2]-1 )
    return false;
  
  const Types::Coordinate from[3] = { idxX*this->m_Delta[0], idxY*this->m_Delta[1], idxZ*this->m_Delta[2] };
  const Types::Coordinate to[3] = { from[0]+this->m_Delta[0], from[1]+this->m_Delta[1], from[2]+this->m_Delta[2] };
  result = this->TrilinearInterpolation( dataPtr, idxX, idxY, idxZ, l, from, to );
  return true;
}

inline bool
UniformVolume::ProbeNoXform
( ProbeInfo& probeInfo, const Self::CoordinateVectorType& location ) const
{
  Self::CoordinateVectorType l( location );
  l -= this->m_Offset;

  if ( (l[0] < 0) || (l[1] < 0) || (l[2] < 0) )
    return false;
  
  const Types::GridIndexType idxX=(Types::GridIndexType) (l[0]/this->m_Delta[0]);
  if ( idxX>=this->m_Dims[0]-1 )
    return false;

  const Types::GridIndexType idxY=(Types::GridIndexType) (l[1]/this->m_Delta[1]);
  if ( idxY>=this->m_Dims[1]-1 )
    return false;

  const Types::GridIndexType idxZ=(Types::GridIndexType) (l[2]/this->m_Delta[2]);
  if ( idxZ>=this->m_Dims[2]-1 )
    return false;

  const Types::Coordinate from[3] = { idxX*this->m_Delta[0], idxY*this->m_Delta[1], idxZ*this->m_Delta[2] };
  const Types::Coordinate to[3] = { from[0]+this->m_Delta[0], from[1]+this->m_Delta[1], from[2]+this->m_Delta[2] };

  return this->GetTrilinear( probeInfo, idxX, idxY, idxZ, l, from, to );
}

inline bool
UniformVolume::FindVoxel
( const Self::CoordinateVectorType& location, Types::GridIndexType *const idx, Types::Coordinate *const from, Types::Coordinate *const to ) const
{
  Self::CoordinateVectorType l( location );
  l -= this->m_Offset;

  if ( (l[0] < 0) || (l[1] < 0) || (l[2] < 0) )
    return false;
  
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<Types::GridIndexType>( l[dim] / this->m_Delta[dim] );
    if ( idx[dim]>=(this->m_Dims[dim]-1) ) 
      return false;
    (to[dim] = (from[dim] = this->m_Offset[dim] + (idx[dim] * this->m_Delta[dim]))) += this->m_Delta[dim];
    }
  
  return true;
}

inline bool
UniformVolume::FindVoxel
( const Self::CoordinateVectorType& location, Types::GridIndexType *const idx ) const
{
  Self::CoordinateVectorType l( location );
  l -= this->m_Offset;
  
  if ( (l[0] < 0) || (l[1] < 0) || (l[2] < 0) )
    return false;
  
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<Types::GridIndexType>( l[dim] / this->m_Delta[dim] );
    if ( idx[dim]>=(this->m_Dims[dim]-1) ) 
      return false;
    }
  return true;
}

inline void
UniformVolume::GetVoxelIndexNoBounds
( const Self::CoordinateVectorType& location, Types::GridIndexType *const idx ) const
{
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<Types::GridIndexType>( floor( (location[dim]-this->m_Offset[dim]) / this->m_Delta[dim]) );
    }
}

inline bool
UniformVolume::FindVoxelByIndex
( const Self::CoordinateVectorType& fracIndex, Types::GridIndexType *const idx, Types::Coordinate *const frac ) const
{
  if ( (fracIndex[0]<0) || (fracIndex[1]<0) || (fracIndex[2]<0) )
    return false;
  
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<Types::GridIndexType>( fracIndex[dim] );
    if ( (idx[dim] >= (this->m_Dims[dim]-1)) ) 
      return false;
    frac[dim] = fracIndex[dim] - idx[dim];
  }
  
  return true;
}

} // namespace cmtk
