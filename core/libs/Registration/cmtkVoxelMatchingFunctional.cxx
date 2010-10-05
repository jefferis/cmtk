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

#include <Registration/cmtkVoxelMatchingFunctional.h>

#include <Base/cmtkVector.h>

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
  
  this->FloatingDims = this->FloatingGrid->GetDims();
  this->FloatingSize = this->FloatingGrid->Size;

  this->m_FloatingCropRegionCoordinates = FloatingGrid->GetHighResCropRegion();
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    this->FloatingInverseDelta[dim] = 1.0 / FloatingGrid->m_Delta[dim];
    this->m_FloatingCropRegionFractional.From()[dim] = this->m_FloatingCropRegionCoordinates.From()[dim] * FloatingInverseDelta[dim];
    this->m_FloatingCropRegionFractional.To()[dim] = this->m_FloatingCropRegionCoordinates.To()[dim] * FloatingInverseDelta[dim];
    }
  
  FloatingDataClass = floating->GetData()->GetDataClass();
}

void
VoxelMatchingFunctional::InitReference( UniformVolume::SmartPtr& reference )
{
  ReferenceGrid = reference;

  this->ReferenceDims = this->ReferenceGrid->GetDims();
  this->ReferenceSize = this->ReferenceGrid->Size;
  this->m_ReferenceCropRegion = ReferenceGrid->CropRegion();

  for ( int dim = 0; dim < 3; ++dim )
    this->ReferenceInvDelta[dim] = 1.0 / ReferenceGrid->m_Delta[dim];

  ReferenceDataClass = reference->GetData()->GetDataClass();
}

const DataGrid::RegionType
VoxelMatchingFunctional::GetReferenceGridRange
( const Vector3D& fromVOI, const Vector3D& toVOI )
{
  DataGrid::IndexType from, to;
  for ( int i = 0; i < 3; ++i )
    {
    from[i] = std::max( this->m_ReferenceCropRegion.From()[i], static_cast<int>( fromVOI[i] * this->ReferenceInvDelta[i] ) );
    to[i] = 1+std::min( this->m_ReferenceCropRegion.To()[i]-1, 1+static_cast<int>( toVOI[i] * this->ReferenceInvDelta[i] ) );
    }

  return DataGrid::RegionType( from, to );
}

} // namespace cmtk
