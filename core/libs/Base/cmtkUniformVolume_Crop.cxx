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

#include "cmtkUniformVolume.h"

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
UniformVolume::SetHighResCropRegion
( const Self::CoordinateRegionType& crop )
{
  if ( !this->m_HighResCropRegion )
    this->m_HighResCropRegion = Self::CoordinateRegionType::SmartPtr( new CoordinateRegionType );

  *this->m_HighResCropRegion = crop;

  for ( int dim = 0; dim<3; ++dim )
    {
    this->CropRegion().From()[dim] = std::max<Self::IndexType::ValueType>( static_cast<Self::IndexType::ValueType>( (crop.From()[dim] - this->m_Offset[dim]) / this->m_Delta[dim] ), 0 );
    this->CropRegion().To()[dim] = 1 + std::min<Self::IndexType::ValueType>( static_cast<Self::IndexType::ValueType>( (crop.To()[dim] - this->m_Offset[dim]) / this->m_Delta[dim] ), this->m_Dims[dim]-1 );
    }
}

const UniformVolume::CoordinateRegionType
UniformVolume::GetHighResCropRegion
() const
{
  if ( this->m_HighResCropRegion )
    {
    return *this->m_HighResCropRegion;
    }
  else
    {
    UniformVolume::CoordinateRegionType region;
    
    for ( int dim = 0; dim<3; ++dim )
      {
      region.From()[dim] = this->m_Offset[dim] + this->m_Delta[dim] * (this->CropRegion().From()[dim]); // take a half pixel off to move between pixels
      region.To()[dim] = this->m_Offset[dim] + this->m_Delta[dim] * (this->CropRegion().To()[dim]-1); // add a hald pixel to move between pixels
      }
    return region;
    }
}

UniformVolume::SmartPtr
UniformVolume::GetCroppedVolume() const
{
  return this->GetCroppedVolume( this->CropRegion() );
}

UniformVolume::SmartPtr
UniformVolume::GetCroppedVolume( const Self::RegionType& region ) const
{
  const Self::IndexType cropDims = region.To() - region.From();
  
  Self::CoordinateVectorType cropSize( cropDims );
  for ( size_t i = 0; i < 3; ++i )
    (cropSize[i] -= 1) *= this->m_Delta[i];
  
  Self::SmartPtr volume( new UniformVolume( cropDims, cropSize ) );
  
  // get cropped data.
  TypedArray::SmartPtr croppedData( this->GetRegionData( region ) );
  volume->SetData( croppedData );

  // prepare new index-to-physical transformation.
  volume->m_IndexToPhysicalMatrix = this->m_IndexToPhysicalMatrix;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      volume->m_IndexToPhysicalMatrix[3][i] += region.From()[j] * volume->m_IndexToPhysicalMatrix[j][i];
  
  // use m_Offset to keep track of new volume origin
  Self::CoordinateVectorType volumeOffset = this->m_Offset;
  for ( int i = 0; i < 3; ++i )
    volumeOffset[i] += (region.From()[i] * this->m_Delta[i]);
  volume->SetOffset( volumeOffset );

  if ( this->m_HighResCropRegion )
    volume->SetHighResCropRegion( *this->m_HighResCropRegion );
  
  volume->CopyMetaInfo( *this, META_IMAGE_ORIENTATION );
  volume->CopyMetaInfo( *this, META_IMAGE_ORIENTATION_ORIGINAL );

  volume->CopyMetaInfo( *this, META_SPACE );

  return volume;
}

} // namespace cmtk
