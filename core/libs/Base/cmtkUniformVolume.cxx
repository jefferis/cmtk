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
#include <cmtkAnatomicalOrientation.h>

#include <cmtkAffineXform.h>
#include <cmtkEigenSystemSymmetricMatrix3x3.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

UniformVolume::UniformVolume
( const DataGrid::IndexType& dims, const Self::CoordinateVectorType& size, TypedArray::SmartPtr& data )
{
  this->SetData( data );
  this->SetDims( dims );
  
  for ( int i=0; i<3; ++i ) {
    Size[i] = static_cast<Types::Coordinate>( size[i] );
    this->m_Delta[i] = ( this->m_Dims[i] == 1 ) ? 0 : Size[i] / (this->m_Dims[i] - 1);
  }

  this->CropRegion() = this->GetWholeImageRegion();
  this->CreateDefaultIndexToPhysicalMatrix();
}

UniformVolume::UniformVolume
( const DataGrid::IndexType& dims, const Types::Coordinate deltaX, const Types::Coordinate deltaY, const Types::Coordinate deltaZ, TypedArray::SmartPtr& data )
{
  this->SetData( data );
  this->SetDims( dims );

  this->m_Delta[0] = deltaX;
  this->m_Delta[1] = deltaY;
  this->m_Delta[2] = deltaZ;

  for ( int i=0; i<3; ++i )
    Size[i] = this->m_Delta[i] * (this->m_Dims[i]-1);

  this->CropRegion() = this->GetWholeImageRegion();
  this->CreateDefaultIndexToPhysicalMatrix();
}

UniformVolume::UniformVolume 
( const UniformVolume& other, const Types::Coordinate resolution, const bool allowUpsampling ) 
{
  Self::IndexType newDims;
  for ( int dim=0; dim<3; ++dim ) 
    {
    Size[dim] = other.Size[dim];
    int new_dims=(int) (Size[dim]/resolution)+1;
    if ( allowUpsampling || (new_dims<=other.m_Dims[dim]) ) 
      {
      newDims[dim]=new_dims;
      this->m_Delta[dim]=Size[dim]/(new_dims-1);
      } 
    else
      {
      if ( other.m_Dims[dim] == 1 ) 
	{
	this->m_Delta[dim] = Size[dim];
	newDims[dim] = 1;
	} 
      else 
	{
	this->m_Delta[dim] = other.m_Delta[dim];
	newDims[dim] = ((int)(Size[dim]/this->m_Delta[dim])) + 1;
	Size[dim] = (newDims[dim]-1) * this->m_Delta[dim];
	}
      }
    }
  
  this->SetDims( newDims );
  TypedArray::SmartPtr resampledData( this->Resample( other ) );
  this->SetData( resampledData	 );
  
  this->m_IndexToPhysicalMatrix = other.m_IndexToPhysicalMatrix;

  this->SetHighResCropRegion( other.GetHighResCropRegion() );
  this->SetOffset( other.m_Offset );
  this->m_MetaInformation = other.m_MetaInformation;
}

UniformVolume*
UniformVolume::Clone( const bool copyData )
{
  if ( copyData ) 
    {
    return this->Clone();
    } 
  else 
    {
    UniformVolume *result = this->CloneGrid();
    result->SetData( this->GetData() );
    
    return result;
    }
}

UniformVolume*
UniformVolume::Clone() const
{
  UniformVolume *result = this->CloneGrid();
  
  if ( this->GetData() )
    {
    TypedArray::SmartPtr clonedData( this->GetData()->Clone() );
    result->SetData( clonedData );
    }
  else
    {
    result->SetData( TypedArray::SmartPtr::Null );
    }
  
  return result;
}

UniformVolume*
UniformVolume::CloneGrid() const
{
  UniformVolume* clone = new UniformVolume( this->m_Dims, Size );
  clone->SetOffset( this->m_Offset );
  clone->m_MetaInformation = this->m_MetaInformation;
  clone->m_IndexToPhysicalMatrix = this->m_IndexToPhysicalMatrix;
  return clone;
}

const UniformVolume::SmartPtr
UniformVolume::GetReoriented( const char* newOrientation ) const
{
  const std::string curOrientation = this->m_MetaInformation[META_IMAGE_ORIENTATION];
  DataGrid::SmartPtr temp( DataGrid::GetReoriented( newOrientation ) );

  AnatomicalOrientation::PermutationMatrix pmatrix( this->m_Dims, curOrientation, newOrientation );
  FixedVector<3,Types::Coordinate> newSize = pmatrix.GetPermutedArray( this->Size );
  
  UniformVolume::SmartPtr result( new UniformVolume( temp->GetDims(), newSize, temp->GetData() ) );
  result->m_Offset = pmatrix.GetPermutedArray( this->m_Offset );
  result->m_IndexToPhysicalMatrix = pmatrix.GetPermutedMatrix( this->m_IndexToPhysicalMatrix );
  result->m_MetaInformation = temp->m_MetaInformation;
  return result;
}

UniformVolume* 
UniformVolume::GetDownsampled( const int downsample, const bool approxIsotropic ) const
{
  if ( approxIsotropic )
    {
    const Types::Coordinate minDelta = std::min<Types::Coordinate>( this->m_Delta[0], std::min<Types::Coordinate>( this->m_Delta[1], this->m_Delta[2] ) );
    const int downsampleByAxis[3] = { std::max<int>( 1, downsample / std::max<int>( 1, static_cast<int>(this->m_Delta[0] / minDelta) ) ),
				      std::max<int>( 1, downsample / std::max<int>( 1, static_cast<int>(this->m_Delta[1] / minDelta) ) ),
				      std::max<int>( 1, downsample / std::max<int>( 1, static_cast<int>(this->m_Delta[2] / minDelta) ) ) };
    return this->GetDownsampled( downsampleByAxis );
    }
  else
    {
    const int downsampleByAxis[3] = { downsample, downsample, downsample };
    return this->GetDownsampled( downsampleByAxis );
    }
}

UniformVolume* 
UniformVolume::GetDownsampled( const int (&downsample)[3] ) const
{
  DataGrid::SmartPtr newDataGrid( this->DataGrid::GetDownsampled( downsample ) );
  TypedArray::SmartPtr newData = newDataGrid->GetData();

  // create downsample grid
  UniformVolume* dsVolume = new UniformVolume( newDataGrid->GetDims(), downsample[0] * this->m_Delta[0], downsample[1] * this->m_Delta[1], downsample[2] * this->m_Delta[2], newData );

  // compute shift of volume origin
  const Vector3D shift( (downsample[0]-1)*this->m_Delta[0]/2, (downsample[1]-1)*this->m_Delta[1]/2, (downsample[2]-1)*this->m_Delta[2]/2 );
  
  // apply shift to origin
  Vector3D offset( this->m_Offset );
  offset += shift;
  dsVolume->SetOffset( offset );
  
  // set crop region while considering new image offset
  dsVolume->SetHighResCropRegion( this->GetHighResCropRegion() );

  dsVolume->m_MetaInformation = this->m_MetaInformation;
  dsVolume->m_IndexToPhysicalMatrix = this->m_IndexToPhysicalMatrix;

  // apply offset shift to index-to-physical matrix
  for ( int axis = 0; axis < 3; ++axis )
    for ( int i = 0; i < 3; ++i )
      {
      dsVolume->m_IndexToPhysicalMatrix[3][i] += (downsample[i]-1) * dsVolume->m_IndexToPhysicalMatrix[axis][i] / 2;
      // also update voxel size
      dsVolume->m_IndexToPhysicalMatrix[axis][i] *= downsample[i];
      }
  
  return dsVolume;
}

UniformVolume* 
UniformVolume::GetInterleavedSubVolume
( const int axis, const int factor, const int idx ) const
{
  Self::IndexType dims;
  Self::CoordinateVectorType size;

  for ( int dim = 0; dim < 3; ++dim )
    {
    dims[dim] = this->m_Dims[dim];
    size[dim] = this->Size[ dim ];
    }
  dims[axis] = this->m_Dims[axis] / factor;
  if ( this->m_Dims[axis] % factor > idx )
    ++dims[axis];
  size[axis] = (dims[axis]-1) * factor * this->m_Delta[axis];
  
  Vector3D offset( 0, 0, 0 );
  offset[axis] = idx * this->m_Delta[axis];
  
  UniformVolume* volume = new UniformVolume( dims, size );
  volume->SetOffset( offset );
  for ( int i = 0; i < dims[axis]; ++i )
    {
    ScalarImage::SmartPtr slice( this->GetOrthoSlice( axis, idx + i * factor ) );
    volume->SetOrthoSlice( axis, i, slice );
    }
  
  volume->m_MetaInformation = this->m_MetaInformation;

  volume->m_IndexToPhysicalMatrix = this->m_IndexToPhysicalMatrix;
  // update coordinate offset according to sub-volume index
  for ( int i = 0; i < 3; ++i )
    volume->m_IndexToPhysicalMatrix[3][i] += idx * volume->m_IndexToPhysicalMatrix[axis][i];
  // scale direction vector along axis
  for ( int i = 0; i < 3; ++i )
    volume->m_IndexToPhysicalMatrix[axis][i] *= factor;
  
  return volume;
}

UniformVolume* 
UniformVolume::GetInterleavedPaddedSubVolume
( const int axis, const int factor, const int idx ) const
{
  int sDims = this->m_Dims[axis] / factor;
  if ( this->m_Dims[axis] % factor > idx )
    ++sDims;

  UniformVolume* volume = new UniformVolume( this->m_Dims, this->Size );
  (volume->CreateDataArray( this->GetData()->GetType() ))->Fill( 0.0 );
  volume->SetOffset( this->m_Offset );
  for ( int i = 0; i < sDims; ++i )
    {
    const size_t sliceIdx = idx + i * factor;
    ScalarImage::SmartPtr slice( this->GetOrthoSlice( axis, sliceIdx ) );
    volume->SetOrthoSlice( axis, sliceIdx, slice );
    }
  
  volume->m_MetaInformation = this->m_MetaInformation;
  volume->m_IndexToPhysicalMatrix = this->m_IndexToPhysicalMatrix;
  return volume;
}

const UniformVolume::RegionType
UniformVolume::GetGridRange
( const Vector3D& fromVOI, const Vector3D& toVOI ) const
{
  Self::IndexType from, to;

  for ( size_t i = 0; i < 3; ++i )
    {
    from[i] = std::max<IndexType::ValueType>( 0, static_cast<IndexType::ValueType>( (fromVOI[i]-this->m_Offset[i]) / this->m_Delta[i] ) );
    to[i] = 1+std::min( this->m_Dims[i]-1, 1+static_cast<IndexType::ValueType>( (toVOI[i]-this->m_Offset[i]) / this->m_Delta[i] ) );
    }

  return UniformVolume::RegionType( from, to );
}

void
UniformVolume::Mirror ( const int axis )
{
  this->DataGrid::ApplyMirrorPlane( axis );
  
  this->CropRegion().From()[ axis ] = this->m_Dims[ axis ] - 1 - this->CropRegion().From()[ axis ];
  this->CropRegion().To()[ axis ] = this->m_Dims[ axis ] - 1 - this->CropRegion().To()[ axis ];
}

ScalarImage*
UniformVolume::GetOrthoSlice
( const int axis, const unsigned int plane ) const
{
  ScalarImage* sliceImage = DataGrid::GetOrthoSlice( axis, plane );
  sliceImage->SetImageSlicePosition( this->GetPlaneCoord( axis, plane ) );

  Vector3D imageOffset( 0, 0, 0 );
  switch ( axis ) 
    {
    case AXIS_X:
      sliceImage->SetPixelSize( this->GetDelta( AXIS_Y, 0 ), this->GetDelta( AXIS_Z, 0 ) );
      imageOffset.Set( this->GetPlaneCoord(AXIS_X, plane), this->GetPlaneCoord(AXIS_Y, 0), this->GetPlaneCoord(AXIS_Z, 0) );
      sliceImage->SetImageDirectionX( Vector3D( 0, 1, 0 ) );
      sliceImage->SetImageDirectionY( Vector3D( 0, 0, 1 ) );
      break;
    case AXIS_Y:
      sliceImage->SetPixelSize( this->GetDelta( AXIS_X, 0 ), this->GetDelta( AXIS_Z, 0 ) );
      imageOffset.Set( this->GetPlaneCoord(AXIS_X, 0), this->GetPlaneCoord(AXIS_Y, plane), this->GetPlaneCoord(AXIS_Z, 0) );
      sliceImage->SetImageDirectionX( Vector3D( 1, 0, 0 ) );
      sliceImage->SetImageDirectionY( Vector3D( 0, 0, 1 ) );
      break;
    case AXIS_Z:
      sliceImage->SetPixelSize( this->GetDelta( AXIS_X, 0 ), this->GetDelta( AXIS_Y, 0 ) );
      imageOffset.Set( this->GetPlaneCoord(AXIS_X, 0), this->GetPlaneCoord(AXIS_Y, 0), this->GetPlaneCoord(AXIS_Z, plane) );
      sliceImage->SetImageDirectionX( Vector3D( 1, 0, 0 ) );
      sliceImage->SetImageDirectionY( Vector3D( 0, 1, 0 ) );
      break;
    }
  
  sliceImage->SetImageOrigin( imageOffset );
  return sliceImage;
}

ScalarImage*
UniformVolume::GetNearestOrthoSlice
( const int axis, const Types::Coordinate location ) const
{
  return this->GetOrthoSlice( axis, this->GetCoordIndex( axis, location ) );
}

ScalarImage*
UniformVolume::GetOrthoSliceInterp
( const int axis, const Types::Coordinate location ) const
{
  unsigned int baseSliceIndex = this->GetCoordIndex( axis, location );

  const Types::Coordinate baseSliceLocation = this->GetPlaneCoord( axis, baseSliceIndex );
  const Types::Coordinate nextSliceLocation = this->GetPlaneCoord( axis, baseSliceIndex+1 ); 

  // only bother to interpolate if we're more than 1% towards the next slice.
  if ( ((location - baseSliceLocation) / (nextSliceLocation-baseSliceLocation )) < 0.01 )
    return this->GetOrthoSlice( axis, baseSliceIndex );

  // only bother to interpolate if we're more than 1% towards the next slice.
  if ( ((nextSliceLocation - location) / (nextSliceLocation-baseSliceLocation )) < 0.01 )
    return this->GetOrthoSlice( axis, baseSliceIndex+1 );

  ScalarImage* image0 = this->GetOrthoSlice( axis, baseSliceIndex );
  ScalarImage* image1 = this->GetOrthoSlice( axis, baseSliceIndex+1 );

  TypedArray::SmartPtr data0 = image0->GetPixelData();
  TypedArray::SmartPtr data1 = image1->GetPixelData();

  Types::Coordinate weight0 = (nextSliceLocation - location) / (nextSliceLocation - baseSliceLocation );

  Types::DataItem value0, value1;
  for ( unsigned int idx = 0; idx < data0->GetDataSize(); ++idx ) 
    {
    if ( data0->Get( value0, idx ) && data1->Get( value1, idx ) ) 
      {
      data0->Set( weight0 * value0 + (1-weight0) * value1, idx );
      } 
    else
      {
      data0->SetPaddingAt( idx );
      }
    }
  delete image1;
  
  image0->SetImageSlicePosition( location );
  image0->SetImageOrigin( weight0 * image0->GetImageOrigin() + (1-weight0) * image1->GetImageOrigin() );
  
  return image0;
}

void
UniformVolume
::GetPrincipalAxes( Matrix3x3<Types::Coordinate>& directions, Vector3D& centerOfMass ) const
{
  centerOfMass = this->GetCenterOfMass();
  const Types::Coordinate xg = centerOfMass[0];
  const Types::Coordinate yg = centerOfMass[1];
  const Types::Coordinate zg = centerOfMass[2];

  Matrix3x3<Types::Coordinate> inertiaMatrix;

  Types::DataItem ixx = 0, iyy = 0, izz = 0, ixy = 0, iyz = 0, izx = 0;
  for ( int k = 0; k < this->m_Dims[2]; ++k )
    {
    const Types::Coordinate Dz = this->GetPlaneCoord( AXIS_Z, k ) - zg;
    const Types::Coordinate Dz2 = Dz * Dz;
    for ( int j = 0; j < this->m_Dims[1]; ++j )
      {
      const Types::Coordinate Dy = this->GetPlaneCoord( AXIS_Y, j ) - yg;
      const Types::Coordinate Dy2 = Dy * Dy;
      for ( int i = 0; i < this->m_Dims[0]; ++i )
	{
        const Types::Coordinate Dx = this->GetPlaneCoord( AXIS_X, i ) - xg;
	const Types::Coordinate Dx2 = Dx * Dx;

	Types::DataItem v;
	if ( this->GetDataAt( v, i, j, k ) )
          {
          ixx += v * ( Dy2 + Dz2 );
          iyy += v * ( Dz2 + Dx2 );
          izz += v * ( Dx2 + Dy2 );
          
          ixy += v * Dx * Dy;
          iyz += v * Dy * Dz;
          izx += v * Dz * Dx;
          }
	}
      }
    }
  inertiaMatrix[0][0] = ixx;
  inertiaMatrix[0][1] = -ixy;
  inertiaMatrix[0][2] = -izx;
  
  inertiaMatrix[1][0] = -ixy;
  inertiaMatrix[1][1] = iyy;
  inertiaMatrix[1][2] = -iyz;
  
  inertiaMatrix[2][0] = -izx;
  inertiaMatrix[2][1] = -iyz;
  inertiaMatrix[2][2] = izz;

  const EigenSystemSymmetricMatrix3x3<Types::Coordinate> eigensystem( inertiaMatrix );
  for ( int n = 0; n < 3; ++n )
    {
    const Vector3D v = eigensystem.GetNthEigenvector( n );
    for ( int i = 0; i < 3; ++i )
      {
      directions[n][i] = v[i];
      }
    }
  
  // correct for negative determinant
  const Types::Coordinate det = directions.Determinant();
  for ( int i = 0; i < 3; i++ )
    {
    directions[2][i] *= det;
    }
   
  for ( int i = 0; i < 3; i++ )
    {
    const Types::Coordinate norm = sqrt( directions[i][0]*directions[i][0] + directions[i][1]*directions[i][1] + directions[i][2]*directions[i][2] );
    for ( int j = 0; j < 3; j++ )
      {
      directions[i][j] /= norm;
      }
    }
}

void
UniformVolume::CreateDefaultIndexToPhysicalMatrix()
{
  this->m_IndexToPhysicalMatrix = AffineXform::MatrixType::IdentityMatrix;
  for ( int axis = 0; axis < 3; ++axis )
    for ( int i = 0; i < 3; ++i )
      this->m_IndexToPhysicalMatrix[axis][i] *= this->m_Delta[axis];
}

void
UniformVolume
::ChangeCoordinateSpace( const std::string& newSpace )
{
  const std::string currentSpace = this->m_MetaInformation[META_SPACE];
  if ( currentSpace == newSpace )
    return; // nothing to do.

  int axesPermutation[3][3];
  AnatomicalOrientation::GetImageToSpaceAxesPermutation( axesPermutation, newSpace.c_str(), currentSpace.c_str() );

  AffineXform::MatrixType newMatrix;
  for ( int i = 0; i < 4; ++i )
    for ( int j = 0; j < 3; ++j )
      for ( int j2 = 0; j2 < 3; ++j2 )
	if ( axesPermutation[j][j2] )
	  newMatrix[i][j] = axesPermutation[j][j2] * this->m_IndexToPhysicalMatrix[i][j2];
  
  this->m_MetaInformation[META_SPACE] = newSpace;
  this->m_IndexToPhysicalMatrix = newMatrix;
}

std::string
UniformVolume
::GetOrientationFromDirections() const
{
  const AffineXform::MatrixType& matrix = this->m_IndexToPhysicalMatrix;
  char orientationString[4] = { 0,0,0,0 };
  AnatomicalOrientation::GetOrientationFromDirections( orientationString, matrix, this->m_MetaInformation[META_SPACE].c_str() );
  return std::string( orientationString );
}

AffineXform::MatrixType
UniformVolume::GetImageToPhysicalMatrix() const
{
  AffineXform::MatrixType matrix = this->m_IndexToPhysicalMatrix;
  for ( int i = 0; i < 3; ++i )
    for ( int j = 0; j < 3; ++j )
      matrix[i][j] /= this->m_Delta[i];

  return matrix;
}

Vector3D
UniformVolume
::IndexToPhysical( const Types::Coordinate i, const Types::Coordinate j, const Types::Coordinate k ) const
{
  const AffineXform::MatrixType& matrix = this->m_IndexToPhysicalMatrix;
  return Vector3D( i * matrix[0][0] + j * matrix[1][0] + k * matrix[2][0] + matrix[3][0],
		   i * matrix[0][1] + j * matrix[1][1] + k * matrix[2][1] + matrix[3][1],
		   i * matrix[0][2] + j * matrix[1][2] + k * matrix[2][2] + matrix[3][2] );  
}

} // namespace cmtk
