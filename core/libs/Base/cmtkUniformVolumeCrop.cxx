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

#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
UniformVolume::SetCropRegionCoordinates
( const UniformVolume::CoordinateRegionType& crop )
{
  for ( int dim = 0; dim<3; ++dim )
    {
    this->CropRegion().From()[dim] = this->GetCoordIndex( dim, std::max<Types::Coordinate>( crop.From()[dim], 0 ) );
    this->CropRegion().To()[dim] = 1 + this->GetCoordIndex( dim, std::min<Types::Coordinate>( crop.To()[dim], this->Size[dim] ) );
    }
}

const UniformVolume::CoordinateRegionType
UniformVolume::GetCropRegionCoordinates
() const
{
  UniformVolume::CoordinateRegionType region;

  for ( int dim = 0; dim<3; ++dim )
    {
    region.From()[dim] = this->m_Offset[dim] + this->m_Delta[dim] * this->CropRegion().From()[dim];
    region.To()[dim] = this->m_Offset[dim] + this->m_Delta[dim] * (this->CropRegion().To()[dim]-1);
    }
  
  return region;
}

UniformVolume* 
UniformVolume::GetCroppedVolume() const
{
  const UniformVolume::IndexType cropDims = this->CropRegion().To() - this->CropRegion().From();
  const Types::Coordinate cropSize[3] = { (cropDims[0]-1) * this->m_Delta[0], (cropDims[1]-1) * this->m_Delta[1], (cropDims[2]-1) * this->m_Delta[2] };
  UniformVolume* volume = new UniformVolume( cropDims, cropSize );

  // get cropped data.
  TypedArray::SmartPtr croppedData( this->GetCroppedData() );
  volume->SetData( croppedData );

  // prepare new index-to-physical transformation.
  volume->m_IndexToPhysicalMatrix = this->m_IndexToPhysicalMatrix;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      volume->m_IndexToPhysicalMatrix[3][i] += this->CropRegion().From()[j] * volume->m_IndexToPhysicalMatrix[j][i];
  
  // use m_Offset to keep track of new volume origin
  volume->SetOffset( this->m_Offset );
  volume->m_Offset += Vector3D( this->GetCropRegionCoordinates().From().begin() );

  volume->m_MetaInformation[META_IMAGE_ORIENTATION]  = this->m_MetaInformation[META_IMAGE_ORIENTATION];
  volume->m_MetaInformation[META_IMAGE_ORIENTATION_ORIGINAL]  = this->m_MetaInformation[META_IMAGE_ORIENTATION_ORIGINAL];

  volume->m_MetaInformation[META_SPACE]  = this->m_MetaInformation[META_SPACE];

  return volume;
}

} // namespace cmtk
