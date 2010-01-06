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

#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

template<class TData> 
inline bool
UniformVolume::ProbeData 
( TData& result, const TData* dataPtr, const Vector3D& location ) const
{
  result=0;
  
  Vector3D l( location );
  l -= this->m_Offset;

  const int idxX=(int) floor(l[0]/this->m_Delta[0]);
  if ( (idxX < 0) || (idxX>=this->m_Dims[0]-1) )
    return false;

  const int idxY=(int) floor(l[1]/this->m_Delta[1]);
  if ( (idxY < 0) || (idxY>=this->m_Dims[1]-1) )
    return false;

  const int idxZ=(int) floor(l[2]/this->m_Delta[2]);
  if ( (idxZ < 0) || (idxZ>=this->m_Dims[2]-1) )
    return false;

  const Types::Coordinate from[3] = { idxX*this->m_Delta[0], idxY*this->m_Delta[1], idxZ*this->m_Delta[2] };
  const Types::Coordinate to[3] = { from[0]+this->m_Delta[0], from[1]+this->m_Delta[1], from[2]+this->m_Delta[2] };
  result = this->TrilinearInterpolation( dataPtr, idxX, idxY, idxZ, l, from, to );
  return true;
}

template<class TData,class TOutputIterator> 
inline bool
UniformVolume::ProbeData 
( const TOutputIterator& result, const std::vector<TData*>& dataPtr, const Vector3D& location ) const
{
  Vector3D l( location );
  l -= this->m_Offset;

  const Types::Coordinate fracX = l[0]/this->m_Delta[0];
  const int idxX = static_cast<int>( floor(fracX) ) ;
  if ( (idxX<0) || (idxX>=this->m_Dims[0]-1) )
    return false;

  const Types::Coordinate fracY = l[1]/this->m_Delta[1];
  const int idxY = static_cast<int>( floor(fracY) ) ;
  if ( (idxY<0) || (idxY>=this->m_Dims[1]-1) )
    return false;

  const Types::Coordinate fracZ = l[2]/this->m_Delta[2];
  const int idxZ = static_cast<int>( floor(fracZ) ) ;
  if ( (idxZ<0) || (idxZ>=this->m_Dims[2]-1) )
    return false;
  
  this->TrilinearInterpolation( result, dataPtr, idxX, idxY, idxZ, fracX - idxX, fracY - idxY, fracZ - idxZ );
  return true;
}

inline bool
UniformVolume::ProbeNoXform
( ProbeInfo& probeInfo, const Vector3D& location ) const
{
  Vector3D l( location );
  l -= this->m_Offset;

  const int idxX=(int) floor(l[0]/this->m_Delta[0]);
  if ( (idxX<0) || (idxX>=this->m_Dims[0]-1) )
    return false;

  const int idxY=(int) floor(l[1]/this->m_Delta[1]);
  if ( (idxY<0) || (idxY>=this->m_Dims[1]-1) )
    return false;

  const int idxZ=(int) floor(l[2]/this->m_Delta[2]);
  if ( (idxZ<0) || (idxZ>=this->m_Dims[2]-1) )
    return false;

  const Types::Coordinate from[3] = { idxX*this->m_Delta[0], idxY*this->m_Delta[1], idxZ*this->m_Delta[2] };
  const Types::Coordinate to[3] = { from[0]+this->m_Delta[0], from[1]+this->m_Delta[1], from[2]+this->m_Delta[2] };

  return this->GetTrilinear( probeInfo, idxX, idxY, idxZ, l, from, to );
}

inline bool
UniformVolume::FindVoxel
( const Vector3D& location, int *const idx, Types::Coordinate *const from, Types::Coordinate *const to ) const
{
  Vector3D l( location );
  l -= this->m_Offset;

  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<int>( floor(l[dim] / this->m_Delta[dim]) );
    if ( (idx[dim]<0) || (idx[dim]>=(this->m_Dims[dim]-1)) ) 
      return false;
    (to[dim] = (from[dim] = this->m_Offset[dim] + (idx[dim] * this->m_Delta[dim]))) += this->m_Delta[dim];
    }
  
  return true;
}

inline bool
UniformVolume::FindVoxel
( const Vector3D& location, int *const idx ) const
{
  Vector3D l( location );
  l -= this->m_Offset;
  
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<int>( floor(l[dim] / this->m_Delta[dim]) );
    if ( (idx[dim]<0) || (idx[dim]>=(this->m_Dims[dim]-1)) ) 
      return false;
    }
  return true;
}

inline void
UniformVolume::GetVoxelIndexNoBounds
( const Vector3D& location, int *const idx ) const
{
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<int>( floor( (location.XYZ[dim]-this->m_Offset.XYZ[dim]) / this->m_Delta[dim]) );
    }
}

inline bool
UniformVolume::FindVoxelByIndex
( const Vector3D& fracIndex, int *const idx, Types::Coordinate *const frac ) const
{
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<int>( floor(fracIndex[dim]) );
    if ( (idx[dim]<0) || (idx[dim] >= (this->m_Dims[dim]-1)) ) 
      return false;
    frac[dim] = fracIndex[dim] - idx[dim];
  }
  
  return true;
}

inline void
UniformVolume::FindVoxelByIndexUnsafe
( const Vector3D& fracIndex, int *const idx, Types::Coordinate *const frac ) const
{
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<int>( floor(fracIndex[dim]) );
    frac[dim] = fracIndex[dim] - idx[dim];
    }
}

inline void
UniformVolume::FindVoxelUnsafe
( const Vector3D& location, int *const idx, Types::Coordinate *const from, Types::Coordinate *const to ) const
{
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    idx[dim] = static_cast<int>( floor((location[dim]-this->m_Offset[dim]) / this->m_Delta[dim]) );
    (to[dim] = from[dim] = this->m_Offset[dim] + idx[dim] * this->m_Delta[dim]) += this->m_Delta[dim];
    }
}

template<class TAccumulator>
ScalarImage*
UniformVolume::ComputeProjection( const int axis ) const
{
  ScalarImage* projectImage = DataGrid::ComputeProjection<TAccumulator>( axis );
  switch ( axis ) 
    {
    case AXIS_X:
      projectImage->SetPixelSize( this->GetDelta( AXIS_Y, 0 ), this->GetDelta( AXIS_Z, 0 ) );
      break;
    case AXIS_Y:
      projectImage->SetPixelSize( this->GetDelta( AXIS_X, 0 ), this->GetDelta( AXIS_Z, 0 ) );
      break;
    case AXIS_Z:
      projectImage->SetPixelSize( this->GetDelta( AXIS_X, 0 ), this->GetDelta( AXIS_Y, 0 ) );
      break;
    }
  return projectImage;
}

} // namespace cmtk
