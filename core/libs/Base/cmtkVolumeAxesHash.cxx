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

#include <cmtkVolumeAxesHash.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

VolumeAxesHash::VolumeAxesHash
( const UniformVolume& volume, const AffineXform* xform, const Types::Coordinate* deltas, const Types::Coordinate* otherOrigin )
{
  // define volume corners
  Vector3D dX(1,0,0), dY(0,1,0), dZ(0,0,1);
  Vector3D V(volume.m_Offset);
  
  dX += volume.m_Offset;
  dY += volume.m_Offset;
  dZ += volume.m_Offset;
    
  if ( xform ) 
    {
    xform->ApplyInPlace(V);
    xform->ApplyInPlace(dX);
    xform->ApplyInPlace(dY);
    xform->ApplyInPlace(dZ);
    }

  dX -= V;
  dY -= V;
  dZ -= V;
  
  if ( otherOrigin )
    {
    V -= Vector3D( otherOrigin );
    }
  
  // Apply post-transformation scaling
  if ( deltas ) 
    {
    Vector3D deltasV( deltas );
    dX = Vector3D::CoordDiv( dX, deltasV );
    dY = Vector3D::CoordDiv( dY, deltasV );
    dZ = Vector3D::CoordDiv( dZ, deltasV );
    V = Vector3D::CoordDiv( V, deltasV );
    }
  
  this->MakeHash( volume, V, dX, dY, dZ );
}

VolumeAxesHash::VolumeAxesHash
( const UniformVolume& volume, const InfinitePlane& mirrorPlane, const Types::Coordinate* deltas )
{
  // define volume corners
  Vector3D dX(1,0,0), dY(0,1,0), dZ(0,0,1);
  Vector3D V(volume.m_Offset);
    
  // apply mirror transformation  
  mirrorPlane.MirrorInPlace(V);
  mirrorPlane.MirrorInPlace(dX);
  dX -= V;
  mirrorPlane.MirrorInPlace(dY);
  dY -= V;
  mirrorPlane.MirrorInPlace(dZ);
  dZ -= V;
  
  // Apply post-transformation scaling
  if ( deltas ) 
    {
    Vector3D deltasV( deltas );
    dX = Vector3D::CoordDiv( dX, deltasV );
    dY = Vector3D::CoordDiv( dY, deltasV );
    dZ = Vector3D::CoordDiv( dZ, deltasV );
    V = Vector3D::CoordDiv( V, deltasV );
    }

  this->MakeHash( volume, V, dX, dY, dZ );
}

void
VolumeAxesHash::MakeHash
( const UniformVolume& volume, const Vector3D& offset, const Vector3D& dX, const Vector3D& dY, const Vector3D& dZ )
{
  int DimsX = volume.m_Dims[0], DimsY = volume.m_Dims[1], DimsZ = volume.m_Dims[2];
  
  this->m_Hash = Memory::AllocateArray<Vector3D*>( 3 );
  assert( this->m_Hash != NULL );

  for ( int dim = 0; dim<3; ++dim ) 
    {
    this->m_Hash[dim] = Memory::AllocateArray<Vector3D>( volume.m_Dims[dim] );
    assert( this->m_Hash[dim] != NULL );
    }

  const Types::Coordinate deltaX = volume.m_Delta[0];
  const Types::Coordinate deltaY = volume.m_Delta[1];
  const Types::Coordinate deltaZ = volume.m_Delta[2];
  
  int idx;
  for ( idx=0; idx<DimsX; ++idx )
    this->m_Hash[0][idx] = deltaX*idx*dX;
  for ( idx=0; idx<DimsY; ++idx )
    this->m_Hash[1][idx] = deltaY*idx*dY;
  for ( idx=0; idx<DimsZ; ++idx )
    (this->m_Hash[2][idx] = deltaZ*idx*dZ) += offset;
}

VolumeAxesHash::~VolumeAxesHash()
{
  assert( this->m_Hash != NULL );

  for ( int dim = 0; dim<3; ++dim ) 
    {
    assert( this->m_Hash[dim] != NULL );
    delete[] this->m_Hash[dim];
    }
  delete[] this->m_Hash;
}

} // namespace cmtk
