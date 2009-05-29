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

#include <cmtkUniformVolume.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
UniformVolume::SetCropRegion
( const Types::Coordinate* cropFromReal, const Types::Coordinate* cropToReal )
{
  if ( cropFromReal ) 
    {
    for ( int dim = 0; dim<3; ++dim )
      {
		  this->CropFromReal[dim] = std::max<Types::Coordinate>( cropFromReal[dim], this->m_Origin[dim] );
      }
    } 
  else
    {
    for ( int dim = 0; dim<3; ++dim )
      {
      this->CropFromReal[dim] = this->m_Origin[dim];
      }
    }

  {
  for ( int dim = 0; dim<3; ++dim )
    this->CropFrom[dim] = this->GetCoordIndex( dim, this->CropFromReal[dim] + this->m_Origin[dim] );
  }

  if ( cropToReal ) 
    {
    for ( int dim = 0; dim<3; ++dim )
      {
      this->CropToReal[dim] = std::min<Types::Coordinate>( cropToReal[dim], this->Size[dim] + this->m_Origin[dim] );
      }
    } 
  else
    {
    for ( int dim = 0; dim<3; ++dim )
      {
      this->CropToReal[dim] = this->Size[dim] + this->m_Origin[dim];
      }
    }
  { 
  for ( int dim = 0; dim<3; ++dim )
    this->CropTo[dim] = 1 + this->GetCoordIndex( dim, this->CropToReal[dim] + this->m_Origin[dim] );
  }
}

void
UniformVolume::SetCropRegion( const int* cropFrom, const int* cropTo )
{
  if ( cropFrom ) 
    {
    for ( int dim = 0; dim < 3; ++dim )
      {
      this->CropFrom[dim] = std::max<int>( 0, std::min<int>( this->Dims[dim], cropFrom[dim] ) );
      }
    } 
  else 
    {
    memset( this->CropFrom, 0, sizeof( this->CropFrom ) );
    }

  { 
  for ( int dim = 0; dim<3; ++dim )
    this->CropFromReal[dim] = this->GetPlaneCoord( dim, this->CropFrom[dim] ) - this->m_Origin[dim];
  }

  if ( cropTo ) 
    {
    for ( int dim = 0; dim < 3; ++dim )
      {
      this->CropTo[dim] = std::max<int>( this->CropFrom[dim], std::min<int>( this->Dims[dim], cropTo[dim] ) );
      }
    } 
  else
    {
    for ( int dim = 0; dim<3; ++dim )
      this->CropTo[dim] = this->Dims[dim];
  }
  { 
  for ( int dim = 0; dim<3; ++dim )
    this->CropToReal[dim] = this->GetPlaneCoord( dim, this->CropTo[dim]-1 ) - this->m_Origin[dim];
  }
}

void
UniformVolume::SetCropRegionFrom( const int* cropFrom )
{
  memcpy( this->CropFrom, cropFrom, sizeof( this->CropFrom ) );

  for ( int dim = 0; dim<3; ++dim )
    this->CropFromReal[dim] = this->GetPlaneCoord( dim, this->CropFrom[dim] ) - this->m_Origin[dim];
}

void
UniformVolume::SetCropRegionTo( const int* cropTo )
{
  memcpy( this->CropTo, cropTo, sizeof( this->CropTo ) );

  for ( int dim = 0; dim<3; ++dim )
    this->CropToReal[dim] = this->GetPlaneCoord( dim, this->CropTo[dim]-1 ) - this->m_Origin[dim];
}

void
UniformVolume::SetCropRegion( const Rect3D* crop )
{
  const int c0[3] = { crop->startX, crop->startY, crop->startZ };
  const int c1[3] = { crop->endX, crop->endY, crop->endZ };

  this->SetCropRegion( c0, c1 );
}

void
UniformVolume::SetCropRegion( const CoordinateRect3D* crop )
{
  const Types::Coordinate c0[3] = { crop->startX, crop->startY, crop->startZ };
  const Types::Coordinate c1[3] = { crop->endX, crop->endY, crop->endZ };

  this->SetCropRegion( c0, c1 );
}

void
UniformVolume::GetCropRegion
( Types::Coordinate *const cropFromReal, Types::Coordinate *const cropToReal ) const
{
  memcpy( cropFromReal, CropFromReal, sizeof( CropFromReal ) );
  memcpy( cropToReal, CropToReal, sizeof( CropToReal ) );
}

void
UniformVolume::GetCropRegion
( int *const cropFrom, int *const cropTo, int *const increments ) const
{
  memcpy( cropFrom, CropFrom, sizeof( CropFrom ) );
  memcpy( cropTo, CropTo, sizeof( CropTo ) );
  
  if ( increments ) 
    {
    increments[0] = cropFrom[0] + Dims[0] * ( cropFrom[1] + Dims[1] * cropFrom[2] );
    increments[1] = cropFrom[0] + (Dims[0] - cropTo[0]);
    increments[2] = Dims[0] * ( cropFrom[1] + ( Dims[1] - cropTo[1]) );
    }
}

int
UniformVolume::GetCropRegionNumVoxels() const 
{
  return (CropTo[0] - CropFrom[0]) * (CropTo[1] - CropFrom[1]) * (CropTo[2] - CropFrom[2]);
}

void
UniformVolume::AutoCrop
( const Types::DataItem threshold, const bool recrop, const int margin )
{
  const TypedArray* data = this->GetData();
  
  int cropFrom[3], cropTo[3];

  if ( ! recrop ) 
    {
    for ( int dim = 0; dim < 3; ++dim ) 
      {
      CropFrom[dim] = 0;
      CropTo[dim] = Dims[dim];
      }
    }
  
  const size_t nextRow = Dims[0] - CropTo[0] + CropFrom[0];
  const size_t nextPlane = Dims[0] * (Dims[1] - CropTo[1] + CropFrom[1]);
  
  bool first = true;
  size_t offset = 
    CropFrom[0] + Dims[0] * ( CropFrom[1] + Dims[1] * CropFrom[2] );

  for ( int z = CropFrom[2]; z < CropTo[2]; ++z, offset += nextPlane )
    for ( int y = CropFrom[1]; y < CropTo[1]; ++y, offset += nextRow )
      for ( int x = CropFrom[0]; x < CropTo[0]; ++x, ++offset ) 
	{
	Types::DataItem value = 0;
	if ( ! data->Get( value, offset ) || (value < threshold) ) continue;
	
	if ( first ) 
	  {
	  cropFrom[0] = cropTo[0] = x;
	  cropFrom[1] = cropTo[1] = y;
	  cropFrom[2] = cropTo[2] = z;
	  first = false;
	  } 
	else
	  {
	  cropFrom[0] = std::min<int>( cropFrom[0], x );
	  cropFrom[1] = std::min<int>( cropFrom[1], y );
	  cropFrom[2] = std::min<int>( cropFrom[2], z );
	  
	  cropTo[0] = std::max<int>( cropTo[0], x );
	  cropTo[1] = std::max<int>( cropTo[1], y );
	  cropTo[2] = std::max<int>( cropTo[2], z );
	  }
	}
  
  for ( int dim = 0; dim < 3; ++dim )
    ++cropTo[dim];
  
  if ( margin ) 
    {
    for ( int dim = 0; dim < 3; ++dim ) 
      {
      cropFrom[dim] = std::max<int>( 0, cropFrom[dim] - margin );
      cropTo[dim] = std::min<int>( Dims[dim], cropTo[dim] + margin );
      }
    }
  
  this->SetCropRegion( cropFrom, cropTo );
}

void
UniformVolume::FillCropBackground( const Types::DataItem value )
{
  const size_t planeSize = Dims[0] * Dims[1];

  size_t offset = CropFrom[2] * planeSize;
  Data->BlockSet( value, 0, offset );

  for ( int z = CropFrom[2]; z < CropTo[2]; ++z ) 
    {
    size_t planeOffset = offset + CropFrom[1] * Dims[0];
    Data->BlockSet( value, offset, planeOffset );

    offset = planeOffset;
    for ( int y = CropFrom[1]; y < CropTo[1]; ++y, offset += Dims[0] ) 
      {
      Data->BlockSet( value, offset, offset+CropFrom[0] );
      Data->BlockSet( value, offset+CropTo[0], offset+Dims[0] );
      }
    
    planeOffset = offset + (Dims[1]-CropTo[1]) * Dims[0];
    Data->BlockSet( value, offset, planeOffset );
    offset = planeOffset;
    }
  
  Data->BlockSet( value, CropTo[2] * planeSize, Dims[2] * planeSize );
}

TypedArray*
UniformVolume::GetCroppedData() const
{
  const TypedArray* srcData = this->GetData();
  if ( ! srcData ) return NULL;

  TypedArray* cropData = 
    srcData->NewTemplateArray( this->GetCropRegionNumVoxels() );
  
  const size_t lineLength = CropTo[0] - CropFrom[0];
  const size_t nextPlane = Dims[0] * (Dims[1] - (CropTo[1] - CropFrom[1]));
  
  size_t toOffset = 0;
  size_t fromOffset = 
    CropFrom[0] + Dims[0] * ( CropFrom[1] + Dims[1] * CropFrom[2] );

  for ( int z = CropFrom[2]; z < CropTo[2]; ++z, fromOffset += nextPlane )
    for ( int y = CropFrom[1]; y < CropTo[1]; ++y, fromOffset += Dims[0] ) 
      {
      srcData->BlockCopy( cropData, toOffset, fromOffset, lineLength );
      toOffset += lineLength;
      }
  
  return cropData;
}

UniformVolume* 
UniformVolume::GetCroppedVolume() const
{
  const int cropDims[3] = { CropTo[0]-CropFrom[0], CropTo[1]-CropFrom[1], CropTo[2]-CropFrom[2] };
  const Types::Coordinate cropSize[3] = { (cropDims[0]-1) * Delta[0], (cropDims[1]-1) * Delta[1], (cropDims[2]-1) * Delta[2] };
  UniformVolume* volume = new UniformVolume( cropDims, cropSize );

  // get cropped data.
  TypedArray::SmartPtr croppedData( this->GetCroppedData() );
  volume->SetData( croppedData );

  // prepare new index-to-physical transformation.
  volume->m_IndexToPhysicalMatrix = this->m_IndexToPhysicalMatrix;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      volume->m_IndexToPhysicalMatrix[3][i] += this->CropFrom[j] * volume->m_IndexToPhysicalMatrix[j][i];
  
  // use m_Origin to keep track of new volume origin
  volume->SetOrigin( this->m_Origin );
  volume->m_Origin += Vector3D( this->CropFromReal );

  volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION]  = this->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION];
  volume->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION_ORIGINAL]  = this->m_MetaInformation[CMTK_META_IMAGE_ORIENTATION_ORIGINAL];

  volume->m_MetaInformation[CMTK_META_SPACE]  = this->m_MetaInformation[CMTK_META_SPACE];

  return volume;
}

} // namespace cmtk
