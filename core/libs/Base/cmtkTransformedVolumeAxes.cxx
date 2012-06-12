/*
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

#include "cmtkTransformedVolumeAxes.h"

#include <cassert>

namespace
cmtk
{

/** \addtogroup Base */
//@{

TransformedVolumeAxes::TransformedVolumeAxes
( const UniformVolume& volume, const AffineXform* xform, const Types::Coordinate* deltas, const Types::Coordinate* otherOrigin )
{
  // define volume corners
  UniformVolume::CoordinateVectorType dX = FixedVectorStaticInitializer<3,Types::Coordinate>::Init(1,0,0);
  UniformVolume::CoordinateVectorType dY = FixedVectorStaticInitializer<3,Types::Coordinate>::Init(0,1,0);
  UniformVolume::CoordinateVectorType dZ = FixedVectorStaticInitializer<3,Types::Coordinate>::Init(0,0,1);
  UniformVolume::CoordinateVectorType V( volume.m_Offset );
  
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
    V -= FixedVector<3,Types::Coordinate>::FromPointer( otherOrigin );
    }
  
  // Apply post-transformation scaling
  if ( deltas ) 
    {
    const UniformVolume::CoordinateVectorType deltasV = UniformVolume::CoordinateVectorType::FromPointer( deltas );
    dX /= deltasV;
    dY /= deltasV;
    dZ /= deltasV;
    V /= deltasV;
    }
  
  this->MakeHash( volume, V, dX, dY, dZ );
}

TransformedVolumeAxes::TransformedVolumeAxes
( const UniformVolume& volume, const ParametricPlane& mirrorPlane, const Types::Coordinate* deltas )
{
  // define volume corners
  UniformVolume::CoordinateVectorType dX = FixedVectorStaticInitializer<3,Types::Coordinate>::Init(1,0,0);
  UniformVolume::CoordinateVectorType dY = FixedVectorStaticInitializer<3,Types::Coordinate>::Init(0,1,0);
  UniformVolume::CoordinateVectorType dZ = FixedVectorStaticInitializer<3,Types::Coordinate>::Init(0,0,1);
  UniformVolume::CoordinateVectorType V( volume.m_Offset );
  
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
    const UniformVolume::CoordinateVectorType deltasV = UniformVolume::CoordinateVectorType::FromPointer( deltas );
    dX /= deltasV;
    dY /= deltasV;
    dZ /= deltasV;
    V /= deltasV;
    }

  this->MakeHash( volume, V, dX, dY, dZ );
}

void
TransformedVolumeAxes::MakeHash
( const UniformVolume& volume, const UniformVolume::SpaceVectorType& offset, const UniformVolume::SpaceVectorType& dX, const UniformVolume::SpaceVectorType& dY, const UniformVolume::SpaceVectorType& dZ )
{
  this->m_Dims = volume.m_Dims;  
  for ( int dim = 0; dim<3; ++dim ) 
    {
    this->m_Hash[dim] = Memory::ArrayC::Allocate<UniformVolume::SpaceVectorType>( this->m_Dims[dim] );
    assert( this->m_Hash[dim] != NULL );
    }

  const Types::Coordinate deltaX = volume.m_Delta[0];
  const Types::Coordinate deltaY = volume.m_Delta[1];
  const Types::Coordinate deltaZ = volume.m_Delta[2];
  
  int idx;
  for ( idx=0; idx < this->m_Dims[0]; ++idx )
    this->m_Hash[0][idx] = deltaX*idx*dX;

  for ( idx=0; idx < this->m_Dims[1]; ++idx )
    this->m_Hash[1][idx] = deltaY*idx*dY;

  for ( idx=0; idx < this->m_Dims[2]; ++idx )
    (this->m_Hash[2][idx] = deltaZ*idx*dZ) += offset;
}

TransformedVolumeAxes::~TransformedVolumeAxes()
{
  for ( int dim = 0; dim<3; ++dim ) 
    {
    assert( this->m_Hash[dim] != NULL );
    Memory::ArrayC::Delete( this->m_Hash[dim] );
    }
}

} // namespace cmtk
