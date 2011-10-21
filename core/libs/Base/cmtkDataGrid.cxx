/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011 SRI International
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

#include "cmtkDataGrid.h"

#include <stdlib.h>

#include <System/cmtkConsole.h>
#include <System/cmtkProgress.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

DataGrid*
DataGrid::CloneVirtual( const bool copyData )
{
  if ( copyData ) 
    {
    return this->CloneVirtual();
    } 
  else 
    {
    DataGrid *result = new Self( this->m_Dims, this->GetData() );
    result->m_CropRegion = this->m_CropRegion;
    
    return result;
    }
}

DataGrid*
DataGrid::CloneVirtual() const
{
  DataGrid *result = new Self( this->m_Dims );
  result->m_CropRegion = this->m_CropRegion;
  
  if ( this->GetData() )
    {
    TypedArray::SmartPtr clonedData( this->GetData()->Clone() );
    result->SetData( clonedData );
    }
  
  return result;
}

const DataGrid::RegionType
DataGrid::GetWholeImageRegion() const
{
  const int zeroes[3] = {0,0,0};
  return Self::RegionType( Self::IndexType( zeroes ), Self::IndexType( this->m_Dims ) );
}

TypedArray::SmartPtr
DataGrid::CreateDataArray( const ScalarDataType dataType, const bool setToZero )
{
  TypedArray::SmartPtr data( TypedArray::Create( dataType, this->GetNumberOfPixels() ) );
  if ( setToZero )
    data->ClearArray();
  this->SetData( data );
  return data;
}

DataGrid* 
DataGrid::GetDownsampledAndAveraged( const int (&downsample)[3] ) const
{
  const int newDims[3] = { (this->m_Dims[0]-1) / downsample[0] + 1, (this->m_Dims[1]-1) / downsample[1] + 1, (this->m_Dims[2]-1) / downsample[2] + 1 };

  DataGrid* newDataGrid = new DataGrid;
  newDataGrid->SetDims( Self::IndexType( newDims ) );
  
  const TypedArray* thisData = this->GetData();
  if ( thisData )
    {
    TypedArray::SmartPtr newData = TypedArray::Create( thisData->GetType(), newDataGrid->GetNumberOfPixels() );

#pragma omp parallel for
    for ( int z = 0; z < newDims[2]; ++z )
      {
      size_t toOffset = z * newDims[0] * newDims[1];
      int oldZ = z * downsample[2];
      int oldY = 0;
      for ( int y = 0; y < newDims[1]; ++y, oldY += downsample[1] ) 
	{
	int oldX = 0;
	for ( int x = 0; x < newDims[0]; ++x, oldX += downsample[0], ++toOffset ) 
	  {
	  Types::DataItem sum = 0;
	  int count = 0;
	  for ( int zz = 0; (zz < downsample[2]) && (oldZ+zz<this->m_Dims[2]); ++zz )
	    {
	    const size_t ofsZ = (oldZ+zz)*this->nextK;
	    for ( int yy = 0; (yy < downsample[1]) && (oldY+yy<this->m_Dims[1]); ++yy )
	      {
	      const size_t ofsYZ = ofsZ + (oldY+yy)*this->nextJ;
	      for ( int xx = 0; (xx < downsample[0]) && (oldX+xx<this->m_Dims[0]); ++xx )
		{
		Types::DataItem value;
		if ( thisData->Get( value, (oldX+xx) + ofsYZ  ) )
		  {
		  sum += value;
		  ++count;
		  }
		}
	      }
	    }
	  if ( count )
	    {
	    newData->Set( sum / count, toOffset );
	    }
	  else
	    {
	    newData->SetPaddingAt( toOffset );
	    }
	  }
	}
      } // end for (z)
    newDataGrid->SetData( TypedArray::SmartPtr( newData ) );
    }

  newDataGrid->CopyMetaInfo( *this, META_IMAGE_ORIENTATION );
  newDataGrid->CopyMetaInfo( *this, META_IMAGE_ORIENTATION_ORIGINAL );
  
  return newDataGrid;
}

DataGrid* 
DataGrid::GetDownsampled( const int (&downsample)[3] ) const
{
  const int newDims[3] = { (this->m_Dims[0]-1) / downsample[0] + 1, (this->m_Dims[1]-1) / downsample[1] + 1, (this->m_Dims[2]-1) / downsample[2] + 1 };

  DataGrid* newDataGrid = new DataGrid;
  newDataGrid->SetDims( Self::IndexType( newDims ) );
  
  const TypedArray* thisData = this->GetData();
  if ( thisData )
    {
    TypedArray::SmartPtr newData = TypedArray::Create( thisData->GetType(), newDataGrid->GetNumberOfPixels() );

#pragma omp parallel for
    for ( int z = 0; z < newDims[2]; ++z )
      {
      size_t toOffset = z * newDims[0] * newDims[1];
      int oldZ = z * downsample[2];
      int oldY = 0;
      for ( int y = 0; y < newDims[1]; ++y, oldY += downsample[1] ) 
	{
	int oldX = 0;
	for ( int x = 0; x < newDims[0]; ++x, oldX += downsample[0], ++toOffset ) 
	  {
	  Types::DataItem value = 0;
	  if ( thisData->Get( value, this->GetOffsetFromIndex( oldX, oldY, oldZ ) ) )
	    newData->Set( value, toOffset );
	  else
	    newData->SetPaddingAt( toOffset );
	  }
	}
      } // end for (z)
    newDataGrid->SetData( TypedArray::SmartPtr( newData ) );
    }

  newDataGrid->CopyMetaInfo( *this, META_IMAGE_ORIENTATION );
  newDataGrid->CopyMetaInfo( *this, META_IMAGE_ORIENTATION_ORIGINAL );
  
  return newDataGrid;
}

const DataGrid::SmartPtr
DataGrid::GetReoriented( const char* newOrientation ) const
{
  std::string curOrientation = this->GetMetaInfo( META_IMAGE_ORIENTATION );
  if ( curOrientation.length() != 3 ) 
    {
    curOrientation = std::string( AnatomicalOrientation::ORIENTATION_STANDARD );
    }

  if ( !newOrientation )
    {
    newOrientation = AnatomicalOrientation::ORIENTATION_STANDARD;
    }

  // 1. get a permutation matrix
  AnatomicalOrientation::PermutationMatrix pmatrix( this->m_Dims, curOrientation, newOrientation );
  
  Self::IndexType newDims = pmatrix.GetPermutedArray( this->m_Dims );
  DataGrid* newDataGrid = new DataGrid( newDims );
  
  const TypedArray* oldData = this->GetData();
  if ( oldData )
    {
    newDataGrid->CreateDataArray( oldData->GetType() );  
    TypedArray* newData = newDataGrid->GetData().GetPtr();
    
    if ( oldData->GetPaddingFlag() ) 
      {
      newData->SetPaddingValue( oldData->GetPaddingValue() );
      }
    
    const char* fromPtr = static_cast<const char*>( oldData->GetDataPtr() );
    char* toPtr = static_cast<char*>( newData->GetDataPtr() );
    
    // 2. loop through the data, applying transformation
    const size_t bytesPerPixel = oldData->GetItemSize();

    int fromPoint[3];
    for ( fromPoint[2] = 0; fromPoint[2] < this->m_Dims[2]; fromPoint[2]++ )
      {
      for ( fromPoint[1] = 0; fromPoint[1] < this->m_Dims[1]; fromPoint[1]++ )
	{
	for ( fromPoint[0] = 0; fromPoint[0] < this->m_Dims[0]; fromPoint[0]++, fromPtr += bytesPerPixel )
	  {
	  memcpy( toPtr + bytesPerPixel * pmatrix.NewOffsetFromPoint( fromPoint ), fromPtr, bytesPerPixel );
	  }
	}
      }
    }
  
  newDataGrid->CopyMetaInfo( *this );
  newDataGrid->SetMetaInfo( META_IMAGE_ORIENTATION, newOrientation );

  return Self::SmartPtr( newDataGrid );
}

void
DataGrid::SetDims( const Self::IndexType& dims )
{
  this->m_Dims = dims;
  this->m_CropRegion = this->GetWholeImageRegion();

  nextI = 1;
  nextJ = nextI * this->m_Dims[0];
  nextK = nextJ * this->m_Dims[1];
  nextIJ = nextI + nextJ;
  nextIK = nextI + nextK;
  nextJK = nextJ + nextK;
  nextIJK = nextI + nextJ + nextK;
}

bool
DataGrid::TrilinearInterpolation
( Types::DataItem& value, const int X, const int Y, const int Z,
  const Self::SpaceVectorType& Location, const Types::Coordinate* from, 
  const Types::Coordinate* to ) const
{
  Types::DataItem corners[8];

  const size_t offset = this->GetOffsetFromIndex( X, Y, Z );

  const TypedArray* data = this->GetData();
  bool data_present = data->Get( corners[0], offset );
  
  if ( X<this->m_Dims[0]-1 ) 
    {
    data_present &= data->Get( corners[1] , offset+nextI );
    
    if ( Y<this->m_Dims[1]-1 ) 
      {
      data_present &= data->Get( corners[3], offset+nextIJ );
      
      if ( Z<this->m_Dims[2]-1 )
	data_present &= data->Get( corners[7], offset+nextIJK );
      else
	return false;
      } 
    else
      {
      return false;
      }
    data_present &= data->Get( corners[5], offset+nextIK );
    } 
  else
    {
    return false;
    }
  
  data_present &= data->Get( corners[2], offset+nextJ );
  data_present &= data->Get( corners[6], offset+nextJK );
  data_present &= data->Get( corners[4], offset+nextK );

  if (data_present) 
    {
    const Types::Coordinate deltaX=1/(to[0]-from[0]), deltaY=1/(to[1]-from[1]), deltaZ=1/(to[2]-from[2]);
    
    const Types::Coordinate revX = deltaX*(Location[0]-from[0]);
    const Types::Coordinate revY = deltaY*(Location[1]-from[1]);
    const Types::Coordinate revZ = deltaZ*(Location[2]-from[2]);
    const Types::Coordinate offsX = 1-revX;
    const Types::Coordinate offsY = 1-revY;
    const Types::Coordinate offsZ = 1-revZ;
    
    value = 
      static_cast<Types::DataItem>( offsZ*(offsY*(offsX*corners[0]+revX*corners[1])+ 
				   revY*(offsX*corners[2]+revX*corners[3]))+
			    revZ*(offsY*(offsX*corners[4]+revX*corners[5])+ 
				  revY*(offsX*corners[6]+revX*corners[7])));
    
    return true;
  }

  return false;
}

TypedArray::SmartPtr
DataGrid::GetDataMirrorPlane( const int axis ) const
{
  TypedArray::SmartPtr result( this->GetData()->Clone() );
  this->MirrorPlaneInPlace( *result, this->m_Dims, axis );
  
  return result;
}

void
DataGrid::ApplyMirrorPlane( const int axis )
{
  this->MirrorPlaneInPlace( *(this->GetData()), this->m_Dims, axis );
}

void
DataGrid::MirrorPlaneInPlace
( TypedArray& data, const Self::IndexType& dims, const int axis )
{
  switch ( axis ) 
    {
    case AXIS_X: 
    {
    int offset = 0;
    for ( int z = 0; z < dims[2]; ++z )
      for ( int y = 0; y < dims[1]; ++y, offset += dims[0] )
	data.BlockReverse( offset, dims[0] );
    }
    break;
    case AXIS_Y: 
    {
    size_t zOffset = 0;
    for ( int z = 0; z < dims[2]; ++z, zOffset += dims[0] * dims[1] ) 
      {
      for ( int y = 0; y < (dims[1] / 2); ++y )
	data.BlockSwap( zOffset + y * dims[0], zOffset + (dims[1] - 1 - y) * dims[0], dims[0] );
      }
    }
    break;
    case AXIS_Z: 
    {
    size_t blockSize = dims[0] * dims[1];
    for ( int z = 0; z < (dims[2] / 2); ++z ) 
      {
      data.BlockSwap( z * blockSize, (dims[2] - 1 - z) * blockSize, blockSize );
      }
    }
    break;
    }
}

ScalarImage*
DataGrid::GetOrthoSlice
( const int axis, const unsigned int plane ) const
{
  unsigned int dims[2], depth, incX, incY, incZ;

  switch ( axis ) 
    {
    case AXIS_X:
      dims[0] = this->m_Dims[1];
      dims[1] = this->m_Dims[2];
      depth = this->m_Dims[0];
      incX = this->m_Dims[0];
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = 1;
      break;
    case AXIS_Y:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[2];
      depth = this->m_Dims[1];
      incX = 1;
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = this->m_Dims[0];
      break;
    case AXIS_Z:
    default:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[1];
      depth = this->m_Dims[2];
      incX = 1;
      incY = this->m_Dims[0];
      incZ = this->m_Dims[0] * this->m_Dims[1];
      break;
    }
  
  const TypedArray* data = this->GetData();
  TypedArray::SmartPtr sliceData = TypedArray::Create( data->GetType(), dims[0] * dims[1] );
  if ( data->GetPaddingFlag() ) 
    {
    sliceData->SetPaddingValue( data->GetPaddingValue() );
    }

  //  if ( (plane < 0) || (plane >= depth) ) {
  if ( (plane >= depth) ) 
    { // need not check < 0 for unsigned int
    sliceData->ClearArray( true /* PaddingData */ );
    } 
  else
    {
    const unsigned int itemSize = data->GetItemSize();
    
    unsigned int sliceOffset = 0;
    unsigned int offset = plane * incZ;
    for ( unsigned int y = 0; y < dims[1]; ++y ) 
      {
      unsigned int offsetY = offset + incY;
      for ( unsigned int x = 0; x < dims[0]; ++x, ++sliceOffset ) 
	{
	unsigned int offsetX = offset + incX;
	
	memcpy( sliceData->GetDataPtr( sliceOffset ), data->GetDataPtr( offset ), itemSize );
	offset = offsetX;
	}
      offset = offsetY;
      }
    }
  
  ScalarImage* sliceImage = new ScalarImage( dims[0], dims[1] );
  sliceImage->SetPixelData( TypedArray::SmartPtr( sliceData ) );
  
  return sliceImage;
}

DataGrid*
DataGrid::ExtractSliceVirtual
( const int axis, const int plane ) const
{
  unsigned int dims[2], incX, incY, incZ;

  switch ( axis ) 
    {
    case AXIS_X:
      dims[0] = this->m_Dims[1];
      dims[1] = this->m_Dims[2];
      incX = this->m_Dims[0];
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = 1;
      break;
    case AXIS_Y:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[2];
      incX = 1;
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = this->m_Dims[0];
      break;
    case AXIS_Z:
    default:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[1];
      incX = 1;
      incY = this->m_Dims[0];
      incZ = this->m_Dims[0] * this->m_Dims[1];
      break;
    }
  
  const TypedArray& data = *(this->GetData());
  TypedArray::SmartPtr sliceData = TypedArray::Create( data.GetType(), dims[0] * dims[1] );
  if ( data.GetPaddingFlag() ) 
    {
    sliceData->SetPaddingValue( data.GetPaddingValue() );
    }

  if ( plane >= this->m_Dims[axis] ) 
    { 
    sliceData->ClearArray( true /* PaddingData */ );
    } 
  else
    {
    const unsigned int itemSize = data.GetItemSize();
    
    unsigned int sliceOffset = 0;
    unsigned int offset = plane * incZ;
    for ( unsigned int y = 0; y < dims[1]; ++y ) 
      {
      unsigned int offsetY = offset + incY;
      for ( unsigned int x = 0; x < dims[0]; ++x, ++sliceOffset ) 
	{
	unsigned int offsetX = offset + incX;
	
	memcpy( sliceData->GetDataPtr( sliceOffset ), data.GetDataPtr( offset ), itemSize );
	offset = offsetX;
	}
      offset = offsetY;
      }
    }
  
  IndexType newDims = this->m_Dims;
  newDims[axis] = 1;

  return new DataGrid( newDims, sliceData );
}

void 
DataGrid::SetOrthoSlice
( const int axis, const unsigned int idx, const ScalarImage* slice )
{
  const TypedArray* sliceData = slice->GetPixelData();
  if ( ! sliceData ) return;

  TypedArray::SmartPtr data = this->GetData();
  if ( ! data )
    {
    data = this->CreateDataArray( sliceData->GetType() );
    }

  unsigned int dims[2], depth, incX, incY, incZ;

  switch ( axis ) 
    {
    case AXIS_X:
      dims[0] = this->m_Dims[1];
      dims[1] = this->m_Dims[2];
      depth = this->m_Dims[0];
      incX = this->m_Dims[0];
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = 1;
      break;
    case AXIS_Y:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[2];
      depth = this->m_Dims[1];
      incX = 1;
      incY = this->m_Dims[0] * this->m_Dims[1];
      incZ = this->m_Dims[0];
      break;
    case AXIS_Z:
    default:
      dims[0] = this->m_Dims[0];
      dims[1] = this->m_Dims[1];
      depth = this->m_Dims[2];
      incX = 1;
      incY = this->m_Dims[0];
      incZ = this->m_Dims[0] * this->m_Dims[1];
      break;
    }
  
  if ( (idx < depth) )
    { // need not check > 0 for unsigned int
    unsigned int sliceOffset = 0;
    unsigned int offset = idx * incZ;
    for ( unsigned int y = 0; y < dims[1]; ++y ) 
      {
      unsigned int offsetY = offset + incY;
      for ( unsigned int x = 0; x < dims[0]; ++x, ++sliceOffset ) 
	{
	unsigned int offsetX = offset + incX;
	
	sliceData->BlockCopy( *data, offset, sliceOffset, 1 );
	offset = offsetX;
	}
      offset = offsetY;
      }
    }
}

FixedVector<3,Types::Coordinate>
DataGrid
::GetCenterOfMassGrid() const
{
  FixedVector<3,Types::Coordinate> com( FixedVector<3,Types::Coordinate>::Init( 0 ) );

  double sumOfSamples = 0;
  size_t ofs = 0;
  for ( int z = 0; z < this->m_Dims[2]; ++z )
    for ( int y = 0; y < this->m_Dims[1]; ++y )
      for ( int x = 0; x < this->m_Dims[0]; ++x, ++ofs )
	{
	Types::DataItem value;
	if ( this->GetDataAt( value, x, y, z ) )
	  {
	  const Types::Coordinate pixelCOM[3] = { value * x, value * y, value * z };
	  com += Self::SpaceVectorType( pixelCOM );
	  sumOfSamples += value;
	  }
	}

  com *= (1.0 / sumOfSamples);

  return com;
}

FixedVector<3,Types::Coordinate>
DataGrid
::GetCenterOfMassGrid(  FixedVector<3,Types::Coordinate>& firstOrderMoment ) const
{
  FixedVector<3,Types::Coordinate> com = this->Self::GetCenterOfMassGrid(); // do not use overloaded function
  firstOrderMoment = FixedVector<3,Types::Coordinate>( FixedVector<3,Types::Coordinate>::Init( 0 ) );

  double sumOfSamples = 0;
  size_t ofs = 0;
  for ( int z = 0; z < this->m_Dims[2]; ++z )
    for ( int y = 0; y < this->m_Dims[1]; ++y )
      for ( int x = 0; x < this->m_Dims[0]; ++x, ++ofs )
	{
	Types::DataItem value;
	if ( this->GetDataAt( value, x, y, z ) )
	  {
	  const Types::Coordinate pixelMoment[3] = { value * fabs(x - com[0]), value * fabs(y - com[1]), value * fabs(z - com[2]) };
	  firstOrderMoment += Self::SpaceVectorType( pixelMoment );
	  sumOfSamples += value;
	  }
	}
  
  firstOrderMoment *= (1.0 / sumOfSamples);

  return com;
}

void 
DataGrid::Print() const
{
  StdErr.printf( "DataGrid at %p\n", this );

  StdErr.printf( "\tthis->m_Dims: [%d,%d,%d])\n", this->m_Dims[0], this->m_Dims[1], this->m_Dims[2] );
  StdErr.printf( "\tData array at %p\n", this );
}

} // namespace cmtk
