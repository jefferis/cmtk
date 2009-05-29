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

#include <cmtkVoxelMatchingFunctional.h>

#include <cmtkVector.h>

#include <assert.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

void
VoxelMatchingFunctional::InitFloating( UniformVolume::SmartPtr& floating )
{
  FloatingGrid = floating;
  
  memcpy( FloatingDims, FloatingGrid->GetDims(), sizeof(FloatingDims) );
  memcpy( FloatingSize, FloatingGrid->Size, sizeof(FloatingSize) );
  FloatingGrid->GetCropRegion( FloatingCropFrom, FloatingCropTo );
  for ( int dim = 0; dim < 3; ++dim ) 
    {
      FloatingInverseDelta.XYZ[dim] = 1.0 / FloatingGrid->Delta[dim];
      FloatingCropFromIndex[dim] = FloatingCropFrom[dim] * FloatingInverseDelta.XYZ[dim];
      FloatingCropToIndex[dim] = FloatingCropTo[dim] * FloatingInverseDelta.XYZ[dim];
    }
  
  FloatingDataClass = floating->GetData()->GetDataClass();
}

void
VoxelMatchingFunctional::InitReference( UniformVolume::SmartPtr& reference )
{
  ReferenceGrid = reference;

  memcpy( ReferenceDims, ReferenceGrid->GetDims(), sizeof(ReferenceDims) );
  memcpy( ReferenceSize, ReferenceGrid->Size, sizeof(ReferenceSize) );
  ReferenceGrid->GetCropRegion( ReferenceCropFrom, ReferenceCropTo );

  for ( int dim = 0; dim < 3; ++dim )
    this->ReferenceInvDelta[dim] = 1.0 / ReferenceGrid->Delta[dim];

  ReferenceDataClass = reference->GetData()->GetDataClass();
}

void
VoxelMatchingFunctional::GetReferenceGridRange
( const Vector3D& fromVOI, const Vector3D& toVOI, Rect3D& voi )
{					      
  voi.startX = std::max( this->ReferenceCropFrom[0], static_cast<int>( fromVOI[0] * this->ReferenceInvDelta[0] ) );
  voi.endX = 1+std::min( this->ReferenceCropTo[0]-1, 1+static_cast<int>( toVOI[0] * this->ReferenceInvDelta[0] ) );
  
  voi.startY = std::max( this->ReferenceCropFrom[1], static_cast<int>( fromVOI[1] * this->ReferenceInvDelta[1] ) );
  voi.endY = 1+std::min( this->ReferenceCropTo[1]-1, 1+static_cast<int>( toVOI[1] * this->ReferenceInvDelta[1] ) );
  
  voi.startZ = std::max( this->ReferenceCropFrom[2], static_cast<int>( fromVOI[2] * this->ReferenceInvDelta[2] ) );
  voi.endZ = 1+std::min( this->ReferenceCropTo[2]-1, 1+static_cast<int>( toVOI[2] * this->ReferenceInvDelta[2] ) );
    
  assert( (voi.startX+1) && (voi.endX<=ReferenceDims[0]) );
  assert( (voi.startY+1) && (voi.endY<=ReferenceDims[1]) );
  assert( (voi.startZ+1) && (voi.endZ<=ReferenceDims[2]) );
}

} // namespace cmtk
