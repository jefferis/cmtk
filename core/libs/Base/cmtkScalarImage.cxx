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

#include "cmtkScalarImage.h"

#include <Base/cmtkMatrix.h>
#include <Base/cmtkCubicSpline.h>
#include <Base/cmtkSurfaceNormal.h>

#include <System/cmtkException.h>

#include <stdlib.h>
#include <cassert>
#include <math.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

ScalarImage::ScalarImage() :
  m_PixelData( NULL ),
  HasROI( false )
{
  this->m_Dims[0] = this->m_Dims[1] = 0;
  this->m_NumberOfFrames = 1;
  this->m_ImageSlicePosition = this->m_ImageTiltAngle = 0;
  this->m_FrameToFrameSpacing = 0;
}

ScalarImage::ScalarImage
( const int dimsx, const int dimsy, const int numberOfFrames ) :
  HasROI( false )
{
  this->m_Dims[0] = dimsx;
  this->m_Dims[1] = dimsy;
  this->m_NumberOfFrames = numberOfFrames;
  this->m_FrameToFrameSpacing = 0;

  this->m_PixelSize[0] = this->m_PixelSize[1] = 1;

  this->m_ImageSlicePosition = this->m_ImageTiltAngle = 0;
}

ScalarImage::ScalarImage
( const ScalarImage& source, 
  const int* roiFrom, const int* roiTo ) :
  HasROI( false )
{
  this->SetDims( source.m_Dims );
  this->SetPixelSize( source.GetPixelSize() );

  this->SetNumberOfFrames( source.GetNumberOfFrames() );
  this->SetFrameToFrameSpacing( source.GetFrameToFrameSpacing() );

  this->SetImageOrigin( source.GetImageOrigin() );
  this->SetImageDirectionX( source.GetImageDirectionX() );
  this->SetImageDirectionY( source.GetImageDirectionY() );
  this->SetImageSlicePosition( source.GetImageSlicePosition() );

  if ( roiFrom && roiTo ) 
    {
    const Self::IndexType::ValueType dims[2] = { roiTo[0] - roiFrom[0], roiTo[1] - roiFrom[1] };
    this->SetDims( Self::IndexType( dims ) );
    
    this->m_ImageOrigin += ( roiFrom[0] * source.GetPixelSize( AXIS_X ) * source.GetImageDirectionX() );
    this->m_ImageOrigin += ( roiFrom[1] * source.GetPixelSize( AXIS_Y ) * source.GetImageDirectionY() );
    
    const TypedArray* sourceData = source.GetPixelData();
    if ( sourceData ) 
      {
      this->CreatePixelData( sourceData->GetType() );
      if ( sourceData->GetPaddingFlag() )
	this->m_PixelData->SetPaddingPtr( sourceData->GetPaddingPtr() );
      
      size_t offset = 0;
      for ( int y = roiFrom[1]; y < roiTo[1]; ++y ) 
	{
	sourceData->ConvertSubArray( this->m_PixelData->GetDataPtr( offset ), this->m_PixelData->GetType(), roiFrom[0] + y * source.m_Dims[AXIS_X], this->m_Dims[AXIS_X] );
	offset += this->m_Dims[AXIS_X];
	}
      }
    } 
  else
    { // if we're not cropping, preserve ROI.
    HasROI = source.HasROI;
    ROI = source.ROI;
    if ( source.GetPixelData() )
      this->SetPixelData( TypedArray::SmartPtr( source.GetPixelData()->Clone() ) );
    }
}

ScalarImage::ScalarImage
( const ScalarImage& source, const Self::RegionType& roi ) :
  HasROI( false )
{
  this->SetDims( roi.To() - roi.From() );
  this->SetPixelSize( source.GetPixelSize() );

  this->SetNumberOfFrames( source.GetNumberOfFrames() );
  this->SetFrameToFrameSpacing( source.GetFrameToFrameSpacing() );

  this->SetImageOrigin( source.GetImageOrigin() );
  this->SetImageDirectionX( source.GetImageDirectionX() );
  this->SetImageDirectionY( source.GetImageDirectionY() );
  this->SetImageSlicePosition( source.GetImageSlicePosition() );

  this->m_ImageOrigin += ( roi.From()[0] * source.GetPixelSize( AXIS_X ) * source.GetImageDirectionX() );
  this->m_ImageOrigin += ( roi.From()[1] * source.GetPixelSize( AXIS_Y ) * source.GetImageDirectionY() );
  
  const TypedArray* sourceData = source.GetPixelData();
  if ( sourceData ) 
    {
    this->CreatePixelData( sourceData->GetType() );
    if ( sourceData->GetPaddingFlag() )
      this->m_PixelData->SetPaddingPtr( sourceData->GetPaddingPtr() );
    
    size_t offset = 0;
    for ( int y = roi.From()[1]; y < roi.To()[1]; ++y ) 
      {
      sourceData->ConvertSubArray( this->m_PixelData->GetDataPtr( offset ), this->m_PixelData->GetType(), roi.From()[0] + y * source.m_Dims[AXIS_X], this->m_Dims[0] );
      offset += this->m_Dims[0];
      }
    }
}

bool 
ScalarImage::GetPixelAt
( Types::DataItem& value, const Types::Coordinate i, const Types::Coordinate j ) const
{
  // check if inside valid pixel index range
  if ( (i < 0) || (i >= this->m_Dims[0]-1) ) return false;
  if ( (j < 0) || (j >= this->m_Dims[1]-1) ) return false;

  // compute index
  const Types::Coordinate I = floor( i );
  const Types::Coordinate J = floor( j );

  size_t ofs = static_cast<size_t>( I + this->m_Dims[0] * J );

  Types::DataItem v00, v01, v10, v11;
  const bool success = 
    this->m_PixelData->Get( v00, ofs ) &&
    this->m_PixelData->Get( v10, ofs + 1 ) &&
    this->m_PixelData->Get( v01, ofs + this->m_Dims[0] ) &&
    this->m_PixelData->Get( v11, ofs + this->m_Dims[0] + 1 );

  const Types::Coordinate ii = i - I;
  const Types::Coordinate jj = j - J;

  // if any of the four lookups hit "Padding Data", return.
  if ( ! success ) return false;

  // compute final value by bilinear interpolation
  value = 
    (1.0 - jj) * ( (1.0 - ii) * v00 + ii * v10 ) + 
    jj * ( (1.0 - ii) * v01 + ii * v11 );
  
  return true;
}

ScalarImage::SpaceVectorType 
ScalarImage::GetImageOrigin( const int frame ) const 
{
  Self::SpaceVectorType origin;
  if ( this->m_NumberOfFrames > 1 ) 
    {
    origin = SurfaceNormal( this->m_ImageDirectionX, this->m_ImageDirectionY ).Get();
    origin *= (frame * this->m_FrameToFrameSpacing) / origin.RootSumOfSquares();
    origin += this->m_ImageOrigin;
    } 
  else
    {
    origin = this->m_ImageOrigin;
    }
  return origin;
}

ScalarImage* 
ScalarImage::Clone() const
{
  ScalarImage *newScalarImage = new ScalarImage( this->m_Dims[0], this->m_Dims[1] );

  newScalarImage->SetPixelSize( this->m_PixelSize );
  newScalarImage->SetImageOrigin( this->m_ImageOrigin );
  newScalarImage->SetImageDirectionX( this->m_ImageDirectionX );
  newScalarImage->SetImageDirectionY( this->m_ImageDirectionY );
  newScalarImage->SetImageSlicePosition( this->m_ImageSlicePosition );

  newScalarImage->SetPixelData( TypedArray::SmartPtr( this->m_PixelData->Clone() ) );
  
  return newScalarImage;
}

ScalarImage* 
ScalarImage::Clone( const bool clonePixelData )
{
  ScalarImage *newScalarImage = new ScalarImage( this->m_Dims[0], this->m_Dims[1] );

  newScalarImage->SetPixelSize( this->m_PixelSize );
  newScalarImage->SetImageOrigin( this->m_ImageOrigin );
  newScalarImage->SetImageDirectionX( this->m_ImageDirectionX );
  newScalarImage->SetImageDirectionY( this->m_ImageDirectionY );
  newScalarImage->SetImageSlicePosition( this->m_ImageSlicePosition );
  
  if ( clonePixelData )
    newScalarImage->SetPixelData( TypedArray::SmartPtr( this->m_PixelData->Clone() ) );
  else
    newScalarImage->SetPixelData( this->m_PixelData );

  return newScalarImage;
}

ScalarImage* 
ScalarImage::Downsample
( const int factorX, int factorY, ScalarImage *const target ) const
{
  if ( ! factorY ) factorY = factorX;

  assert( this->m_NumberOfFrames == 1 );

  ScalarImage *newScalarImage = target;
  if ( ! newScalarImage )
    newScalarImage = new ScalarImage( this->m_Dims[0] / factorX, this->m_Dims[1] / factorY );
  
  newScalarImage->SetPixelSize( this->m_PixelSize[0] * factorX, this->m_PixelSize[1] * factorY );
  
  Self::SpaceVectorType imageOrigin( this->m_ImageOrigin );
  imageOrigin += (0.5 * this->m_PixelSize[0] / this->m_ImageDirectionX.RootSumOfSquares()) * this->m_ImageDirectionX;
  imageOrigin += (0.5 * this->m_PixelSize[1] / this->m_ImageDirectionY.RootSumOfSquares()) * this->m_ImageDirectionY;

  newScalarImage->SetImageOrigin( imageOrigin );
  newScalarImage->SetImageDirectionX( this->m_ImageDirectionX );
  newScalarImage->SetImageDirectionY( this->m_ImageDirectionY );
  newScalarImage->SetImageSlicePosition( this->m_ImageSlicePosition );
  
  newScalarImage->CreatePixelData( this->m_PixelData->GetType() );

  // compensate for dimensions that connot be evenly divided by downsampling 
  // factor.
  const int dimsY = (this->m_Dims[1] / factorY) * factorY;
  const int dimsX = (this->m_Dims[0] / factorX) * factorX;
  const Types::DataItem factorXY = 1.0 / (factorX * factorY);

  int j = 0;
  for ( int y = 0; y < dimsY; y += factorY, ++j ) 
    {
    int i = 0;
    for ( int x = 0; x < dimsX; x += factorX, ++i ) 
      {
      Types::DataItem pixel = 0;
      for ( int yy = 0; yy < factorY; ++yy )
	for ( int xx = 0; xx < factorX; ++xx )
	  pixel += this->GetPixelAt( x + xx, y + yy );
      
      newScalarImage->SetPixelAt( i, j, pixel * factorXY );
      }
    }
  return newScalarImage;  
}

TypedArray::SmartPtr
ScalarImage::GetSobelFiltered( const bool horizontal, const bool absolute ) 
  const
{
  if ( !this->m_PixelData )
    throw( Exception( "No image data in ScalarImage::GetSobelFiltered()" ) );

  TypedArray::SmartPtr result = ( absolute ) ?
    TypedArray::Create( GetUnsignedDataType( this->m_PixelData->GetType() ), this->m_PixelData->GetDataSize() ) :
    TypedArray::Create( GetSignedDataType( this->m_PixelData->GetType() ), this->m_PixelData->GetDataSize() );
  
  Types::DataItem fov[3][3];
  size_t offset = 0;

  if ( horizontal )
    for ( int y = 0; y < this->m_Dims[1]; ++y )
      for ( int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
	{
	Types::DataItem value = 0;
	if ( x && y && (x<this->m_Dims[0]-1) && (y<this->m_Dims[1]-1) ) 
	  {
	  for ( byte dy=0; dy<3; ++dy )
	    for ( byte dx=0; dx<3; ++dx )
	      this->m_PixelData->Get( fov[dx][dy], x-1+dx+ this->m_Dims[0] * ( y-1+dy ) );
	  
	  value = fov[2][0] - fov[0][0] + 2 * ( fov[2][1] - fov[0][1] ) + fov[2][2] - fov[0][2];
	  
	  if ( absolute ) value = fabs( value );
	  }
	result->Set( value, offset );
	}
  else
    for ( int y = 0; y < this->m_Dims[1]; ++y )
      for ( int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
	{
	Types::DataItem value = 0;
	if ( x && y && (x<this->m_Dims[0]-1) && (y<this->m_Dims[1]-1) ) 
	  {
	  for ( byte dy=0; dy<3; ++dy )
	    for ( byte dx=0; dx<3; ++dx )
	      this->m_PixelData->Get( fov[dx][dy], x-1+dx+ this->m_Dims[0] * ( y-1+dy ) );
	  
	  value = fov[0][0] - fov[0][2] + 2 * ( fov[1][0] - fov[1][2] ) + fov[2][0] - fov[2][2];
	  
	  if ( absolute ) value = fabs( value );
	  }
	result->Set( value, offset );
	}
  
  return result;
}

ScalarImage* operator- 
( const ScalarImage& image0, const ScalarImage& image1 )
{
  ScalarImage *result = new ScalarImage( image0.m_Dims[0], image0.m_Dims[1] );

  const TypedArray *data0 = image0.GetPixelData();
  const TypedArray *data1 = image1.GetPixelData();

  size_t numberOfPixels = image0.GetNumberOfPixels();

  TypedArray::SmartPtr pixelData( TypedArray::Create( GetSignedDataType( data0->GetType() ), numberOfPixels ) );

  Types::DataItem pixel0, pixel1;
  for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
    {
    if ( data0->Get( pixel0, idx ) && data1->Get( pixel1, idx ) ) 
      {
      pixelData->Set( pixel0 - pixel1, idx );
      } 
    else
      {
      pixelData->SetPaddingAt( idx );
      }
    }
  
  result->SetPixelData( pixelData );
  
  return result;
}

ScalarImage&
ScalarImage::operator-=
( const ScalarImage& source )
{
  TypedArray *data0 = this->GetPixelData().GetPtr();
  const TypedArray *data1 = source.GetPixelData();

  size_t numberOfPixels = this->GetNumberOfPixels();
  Types::DataItem pixel0, pixel1;
  for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
    {
    if ( data0->Get( pixel0, idx ) && data1->Get( pixel1, idx ) ) 
      {
      data0->Set( pixel0 - pixel1, idx );
      } 
    else
      {
      data0->SetPaddingAt( idx );
      }
    }
  
  return *this;
}

void ScalarImage::Mirror( const bool horizontal, const bool vertical )
{
  if ( vertical ) 
    {
    for ( int y = 0; y < this->m_Dims[1]/2; ++y ) 
      {
      this->m_PixelData->BlockSwap( y * this->m_Dims[0], (this->m_Dims[1]-y-1) * this->m_Dims[0], this->m_Dims[0] );
      }
    this->m_ImageOrigin = this->m_ImageOrigin + ((this->m_Dims[1]-1) * this->m_PixelSize[1] / this->m_ImageDirectionY.RootSumOfSquares()) * this->m_ImageDirectionY;
    this->m_ImageDirectionY *= (-1.0);
    }
  
  if ( horizontal ) 
    {
    for ( int y = 0; y < this->m_Dims[1]; ++y ) 
      {
      this->m_PixelData->BlockReverse( y * this->m_Dims[0], this->m_Dims[0] );
      }
    this->m_ImageOrigin = this->m_ImageOrigin + ((this->m_Dims[1]-1) * this->m_PixelSize[0] / this->m_ImageDirectionX.RootSumOfSquares()) * this->m_ImageDirectionX;
    this->m_ImageDirectionX *= (-1.0);
    }
}

void
ScalarImage::AdjustToIsotropic
( const Types::Coordinate pixelSize, const bool interpolate )
{
  if ( pixelSize < this->m_PixelSize[0] )
    {
    // fake pixel size Y, then simply adjust aspect ratio
    const Types::Coordinate savePixelSizeY = this->m_PixelSize[1];
    this->m_PixelSize[1] = pixelSize;
    this->AdjustAspectX( interpolate );
    this->m_PixelSize[1] = savePixelSizeY;    
    }

  // now we can simply adjust aspect again in the other dimension
  if ( this->m_PixelSize[0] < this->m_PixelSize[1] )
    {
    this->AdjustAspectY( interpolate );
    }
}

void
ScalarImage::AdjustAspect( const bool interpolate )
{
  if ( this->m_PixelSize[0] < this->m_PixelSize[1] )
    this->AdjustAspectY( interpolate );
  else
    if ( this->m_PixelSize[0] > this->m_PixelSize[1] )
      this->AdjustAspectX( interpolate );
}

void
ScalarImage::AdjustAspectY( const bool interpolate )
{
  if ( this->m_Dims[0] < 2 )
    return;
  
  const int newDimsY = static_cast<int>( (this->m_Dims[1]-1) * this->m_PixelSize[1]/this->m_PixelSize[0] ) + 1;
  
  TypedArray::SmartPtr scaled = TypedArray::Create( this->m_PixelData->GetType(), this->m_Dims[0] * newDimsY );
  
  if ( interpolate ) 
    {
    // with interpolation
    std::vector<Types::DataItem> row0( this->m_Dims[0] ), row1( this->m_Dims[0] );
    this->m_PixelData->GetSubArray( &(row0[0]), 0, this->m_Dims[0] );
    this->m_PixelData->GetSubArray( &(row1[0]), this->m_Dims[0], this->m_Dims[0] );
    
    Types::Coordinate scanLine = 0;
    int ySource = 0;
    size_t offset = 0;
    for ( int y = 0; y < newDimsY; ++y ) 
      {
      Types::Coordinate factor = scanLine / this->m_PixelSize[1];
      for ( int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
	{
	scaled->Set( (1.0 - factor ) * row0[x] + factor * row1[x], offset );
	}
      
      scanLine += this->m_PixelSize[0];
      while ( (ySource < this->m_Dims[1]) && (scanLine >= this->m_PixelSize[1]) ) 
	{
	++ySource;
	row0 = row1;
	this->m_PixelData->GetSubArray( &(row1[0]), (1+ySource) * this->m_Dims[0], this->m_Dims[0] );
	scanLine -= this->m_PixelSize[1];
	}
      }
    } 
  else
    {
    // no interpolation; can do with simply block copying
    char *scaledPtr = static_cast<char*>( scaled->GetDataPtr() );
    const char *pixelPtr = static_cast<char*>( this->m_PixelData->GetDataPtr() );
    
    Types::Coordinate scanLineOffset = this->m_PixelSize[1] / 2;
    Types::Coordinate scanLine = scanLineOffset; // correct offset for NN
    int ySource = 0;
    for ( int y = 0; y < newDimsY; ++y ) 
      {
      memcpy( scaledPtr, pixelPtr, scaled->GetItemSize() * this->m_Dims[0] );
      scanLine += this->m_PixelSize[0];
      while ( (ySource < this->m_Dims[1]) && (scanLine >= this->m_PixelSize[1]) ) 
	{
	++ySource;
	pixelPtr += this->m_PixelData->GetItemSize() * this->m_Dims[0];
	scanLine -= this->m_PixelSize[1];
	}
      scaledPtr += scaled->GetItemSize() * this->m_Dims[0];
      }
    }
  
  this->m_PixelSize[1] = this->m_PixelSize[0];
  this->m_Dims[1] = newDimsY;
  this->SetPixelData( scaled );
}

void
ScalarImage::AdjustAspectX( const bool interpolate )
{
  if ( this->m_Dims[1] < 2 )
    return;
  
  const int newDimsX = static_cast<int>( (this->m_Dims[0]-1) * this->m_PixelSize[0]/this->m_PixelSize[1] ) + 1;
  
  TypedArray::SmartPtr scaled = TypedArray::Create( this->m_PixelData->GetType(), newDimsX * this->m_Dims[1] );
  
  if ( interpolate ) 
    {
    std::vector<Types::Coordinate> factor( newDimsX );
    std::vector<int> fromPixel( newDimsX );
    
    Types::Coordinate scanLine = 0;
    int xSource = 0;
    for ( int x = 0; x < newDimsX; ++x ) 
      {
      fromPixel[x] = xSource;
      factor[x] = scanLine / this->m_PixelSize[0];
      scanLine += this->m_PixelSize[1];
      while ( (xSource < this->m_Dims[0]) && (scanLine >= this->m_PixelSize[0]) ) 
	{
	++xSource;
	scanLine -= this->m_PixelSize[0];
	}
      }
    
    std::vector<Types::DataItem> rowFrom( this->m_Dims[0] );
    size_t offset = 0;
    for ( int y = 0; y < this->m_Dims[1]; ++y ) 
      {
      this->m_PixelData->GetSubArray( &(rowFrom[0]), y * this->m_Dims[0], this->m_Dims[0] );
      for ( int x = 0; x < newDimsX; ++x, ++offset ) 
	{
	scaled->Set( (1.0 - factor[x] ) * rowFrom[fromPixel[x]] + factor[x] * rowFrom[fromPixel[x]+1], offset );
	}
      }
    } 
  else
    {
    Types::Coordinate scanLine = this->m_PixelSize[0] / 2; // correct offset for NN
    int xSource = 0;
    std::vector<int> fromPixel( newDimsX );
    for ( int x = 0; x < newDimsX; ++x ) 
      {
      fromPixel[x] = xSource * scaled->GetItemSize();
      scanLine += this->m_PixelSize[1];
      while ( (xSource < this->m_Dims[0]) && (scanLine >= this->m_PixelSize[0]) ) 
	{
	++xSource;
	scanLine -= this->m_PixelSize[0];
	}
      }
    
    // no interpolation; can do with simply block copying
    char *scaledPtr = static_cast<char*>( scaled->GetDataPtr() );
    const char *pixelPtr = static_cast<char*>( this->m_PixelData->GetDataPtr() );
    
    for ( int y = 0; y < this->m_Dims[1]; ++y ) 
      {
      for ( int x = 0; x < newDimsX; ++x ) 
	{
	memcpy( scaledPtr, pixelPtr+fromPixel[x], scaled->GetItemSize() );
	scaledPtr += scaled->GetItemSize();
	}
      pixelPtr += scaled->GetItemSize() * this->m_Dims[0];
      }
    }
  
  this->m_PixelSize[0] = this->m_PixelSize[1];
  this->m_Dims[0] = newDimsX;
  this->SetPixelData( scaled );
}

void
ScalarImage::ProjectPixel
( const Self::SpaceVectorType& v, int& i, int& j ) const
{
  Self::SpaceVectorType p(v);
  p -= this->m_ImageOrigin;
  
  i = MathUtil::Round( ( p * this->m_ImageDirectionX ) / ( this->m_ImageDirectionX.SumOfSquares() * this->m_PixelSize[0] ) );
  j = MathUtil::Round( ( p * this->m_ImageDirectionY ) / ( this->m_ImageDirectionY.SumOfSquares() * this->m_PixelSize[1] ) );
}

void 
ScalarImage::Print() const
{
  StdErr.printf( "ScalarImage at %p\n", this );

  StdErr.indent();
  StdErr.printf( "Dimensions: [%d,%d]\n", this->m_Dims[0], this->m_Dims[1] );
  StdErr.printf( "Pixel size: [%f,%f]\n", this->m_PixelSize[0], this->m_PixelSize[1] );
  StdErr.printf( "Origin: [%f,%f,%f]\n", this->m_ImageOrigin[0], this->m_ImageOrigin[1], this->m_ImageOrigin[2] );

  StdErr.printf( "DirectionX: [%f,%f,%f]\n", this->m_ImageDirectionX[0], this->m_ImageDirectionX[1], this->m_ImageDirectionX[2] );
  StdErr.printf( "DirectionY: [%f,%f,%f]\n", this->m_ImageDirectionY[0], this->m_ImageDirectionY[1], this->m_ImageDirectionY[2] );
  StdErr.unindent();
}

} // namespace cmtk
