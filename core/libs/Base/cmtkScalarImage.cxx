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

#include <cmtkScalarImage.h>

#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include <cmtkMatrix.h>
#include <cmtkCubicSpline.h>

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
( const unsigned int dimsx, const unsigned int dimsy, const unsigned int numberOfFrames ) :
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
( const ScalarImage* other, 
  const unsigned int* roiFrom, const unsigned int* roiTo ) :
  HasROI( false )
{
  this->SetDims( other->GetDims( AXIS_X ), other->GetDims( AXIS_Y ) );
  this->SetPixelSize( other->GetPixelSize() );

  this->SetNumberOfFrames( other->GetNumberOfFrames() );
  this->SetFrameToFrameSpacing( other->GetFrameToFrameSpacing() );

  this->SetImageOrigin( other->GetImageOrigin() );
  this->SetImageDirectionX( other->GetImageDirectionX() );
  this->SetImageDirectionY( other->GetImageDirectionY() );
  this->SetImageSlicePosition( other->GetImageSlicePosition() );

  if ( roiFrom && roiTo ) 
    {
    this->SetDims( roiTo[0] - roiFrom[0], roiTo[1] - roiFrom[1] );
    
    this->m_ImageOrigin += ( roiFrom[0] * other->GetPixelSize( AXIS_X ) * other->GetImageDirectionX() );
    this->m_ImageOrigin += ( roiFrom[1] * other->GetPixelSize( AXIS_Y ) * other->GetImageDirectionY() );
    
    const TypedArray* otherData = other->GetPixelData();
    if ( otherData ) 
      {
      this->CreatePixelData( otherData->GetType() );
      if ( otherData->GetPaddingFlag() )
	this->m_PixelData->SetPaddingPtr( otherData->GetPaddingPtr() );
      
      size_t offset = 0;
      for ( unsigned int y = roiFrom[1]; y < roiTo[1]; ++y ) 
	{
	otherData->ConvertSubArray( this->m_PixelData->GetDataPtr( offset ), this->m_PixelData->GetType(), roiFrom[0] + y * other->GetDims( AXIS_X ), this->m_Dims[0] );
	offset += this->m_Dims[0];
	}
      }
    } 
  else
    { // if we're not cropping, preserve ROI.
    HasROI = other->HasROI;
    ROI = other->ROI;
    if ( other->GetPixelData() )
      this->SetPixelData( TypedArray::SmartPtr( other->GetPixelData()->Clone() ) );
    }
}

ScalarImage::ScalarImage
( const ScalarImage* other, const IntROI2D* roi ) :
  HasROI( false )
{
  this->SetDims( roi->To[0] - roi->From[0], roi->To[1] - roi->From[1] );
  this->SetPixelSize( other->GetPixelSize() );

  this->SetNumberOfFrames( other->GetNumberOfFrames() );
  this->SetFrameToFrameSpacing( other->GetFrameToFrameSpacing() );

  this->SetImageOrigin( other->GetImageOrigin() );
  this->SetImageDirectionX( other->GetImageDirectionX() );
  this->SetImageDirectionY( other->GetImageDirectionY() );
  this->SetImageSlicePosition( other->GetImageSlicePosition() );

  this->m_ImageOrigin += ( roi->From[0] * other->GetPixelSize( AXIS_X ) * other->GetImageDirectionX() );
  this->m_ImageOrigin += ( roi->From[1] * other->GetPixelSize( AXIS_Y ) * other->GetImageDirectionY() );
  
  const TypedArray* otherData = other->GetPixelData();
  if ( otherData ) 
    {
    this->CreatePixelData( otherData->GetType() );
    if ( otherData->GetPaddingFlag() )
      this->m_PixelData->SetPaddingPtr( otherData->GetPaddingPtr() );
    
    size_t offset = 0;
    for ( int y = roi->From[1]; y < roi->To[1]; ++y ) 
      {
      otherData->ConvertSubArray( this->m_PixelData->GetDataPtr( offset ), this->m_PixelData->GetType(), roi->From[0] + y * other->GetDims( AXIS_X ), this->m_Dims[0] );
      offset += this->m_Dims[0];
      }
    }
}

ScalarImage*
ScalarImage::InterpolateFrom
( const ScalarImage* grid, const CoordinateMatrix3x3* matrix, const cmtk::Interpolators::InterpolationEnum interpolation ) const
{
  int dimsX = grid->m_Dims[0], dimsY = grid->m_Dims[1];
  ScalarImage* result = new ScalarImage( dimsX, dimsY );

  // if we don;t have pixel data, don't interpolate
  if ( !this->GetPixelData() ) return result;

  // create pixel data of same type as ours
  result->CreatePixelData( this->GetPixelData()->GetType() );
  TypedArray* data = result->GetPixelData().GetPtr();

  // compute transformed region origin and pixel directions
  Types::Coordinate origin[2] = { 0, 0 };
  Types::Coordinate deltax[2] = { grid->m_PixelSize[0], 0 };
  Types::Coordinate deltay[2] = { 0, grid->m_PixelSize[1] };

  matrix->Multiply( origin );
  origin[0] /= this->m_PixelSize[0];
  origin[1] /= this->m_PixelSize[1];

  matrix->Multiply( deltax );
  deltax[0] /= this->m_PixelSize[0];
  deltax[1] /= this->m_PixelSize[1];
  deltax[0] -= origin[0];
  deltax[1] -= origin[1];

  matrix->Multiply( deltay );
  deltay[0] /= this->m_PixelSize[0];
  deltay[1] /= this->m_PixelSize[1];
  deltay[0] -= origin[0];
  deltay[1] -= origin[1];
  
  Types::DataItem value;
  size_t offset = 0;

  switch ( interpolation )
    {
    default:
    case cmtk::Interpolators::LINEAR:
      for ( int y = 0; y < dimsY; ++y ) 
	{
	Types::Coordinate row[2] = { origin[0], origin[1] };
	for ( int x = 0; x < dimsX; ++x, ++offset ) 
	  {
	  if ( this->GetPixelAt( value, row[0], row[1] ) )
	    data->Set( value, offset );
	  else
	    data->SetPaddingAt( offset );
	  row[0] += deltax[0];
	  row[1] += deltax[1];
	  }
	origin[0] += deltay[0];
	origin[1] += deltay[1];
	}
      break;
    case cmtk::Interpolators::CUBIC:
      for ( int y = 0; y < dimsY; ++y ) 
	{
	Types::Coordinate row[2] = { origin[0], origin[1] };
	for ( int x = 0; x < dimsX; ++x, ++offset ) 
	  {
	  if ( this->GetPixelAtCubic( value, row[0], row[1] ) )
	    data->Set( value, offset );
	  else
	    data->SetPaddingAt( offset );
	  row[0] += deltax[0];
	  row[1] += deltax[1];
	  }
	origin[0] += deltay[0];
	origin[1] += deltay[1];
	}
      break;
    }

  return result;
}

void
ScalarImage::ApplyBinaryMask
( const ScalarImage* maskImage, const bool invert )
{
  const TypedArray* maskData = maskImage->GetPixelData();
  if ( ! maskData ) return;
  if ( ! this->m_PixelData ) return;

  const size_t numberOfPixels = this->m_PixelData->GetDataSize();
  for ( size_t idx = 0; idx < numberOfPixels; ++idx ) 
    {
    Types::DataItem maskValue;
    if ( maskData->Get( maskValue, idx ) ) 
      {
      if ( (maskValue == 0) ^ invert )
	this->m_PixelData->SetPaddingAt( idx );
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

bool 
ScalarImage::GetPixelAtCubic
( Types::DataItem& value, const Types::Coordinate i, const Types::Coordinate j ) const
{
  // check if inside valid pixel index range
  if ( (i < 1) || (i >= this->m_Dims[0]-2) ) return false;
  if ( (j < 1) || (j >= this->m_Dims[1]-2) ) return false;

  // compute index
  const Types::Coordinate I = floor( i );
  const Types::Coordinate J = floor( j );

  size_t ofs = static_cast<size_t>( (I-1) + this->m_Dims[0] * (J-1) );

  const Types::Coordinate ii = (i - I);
  const Types::Coordinate jj = (j - J);

  const Types::Coordinate spI[4] = 
    { CubicSpline::InterpSpline( 0, ii ),
      CubicSpline::InterpSpline( 1, ii ),
      CubicSpline::InterpSpline( 2, ii ),
      CubicSpline::InterpSpline( 3, ii ) };

  const Types::Coordinate spJ[4] = 
    { CubicSpline::InterpSpline( 0, jj ),
      CubicSpline::InterpSpline( 1, jj ),
      CubicSpline::InterpSpline( 2, jj ),
      CubicSpline::InterpSpline( 3, jj ) };

  const TypedArray* data = this->GetPixelData();

  value = 0;
  Types::DataItem item;
  for ( int j=0; j<4; ++j )
    {
    Types::DataItem value_j = 0;
    for ( int i=0; i<4; ++i ) 
      {
      if ( data->Get( item, ofs + i + j * this->m_Dims[0] ) )
	value_j += static_cast<Types::DataItem>( item * spI[i] );
      else
	return false;
      }
    value += static_cast<Types::DataItem>( value_j * spJ[j] );
    }
  
  return true;
}

Vector3D 
ScalarImage::GetImageOrigin( const unsigned int frame ) const 
{
  Vector3D origin;
  if ( this->m_NumberOfFrames > 1 ) 
    {
    origin.SetNormal( this->m_ImageDirectionX, this->m_ImageDirectionY );
    origin *= (frame * this->m_FrameToFrameSpacing) / origin.EuclidNorm();
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
  
  Vector3D imageOrigin( this->m_ImageOrigin );
  imageOrigin += (0.5 * this->m_PixelSize[0] / this->m_ImageDirectionX.EuclidNorm()) * this->m_ImageDirectionX;
  imageOrigin += (0.5 * this->m_PixelSize[1] / this->m_ImageDirectionY.EuclidNorm()) * this->m_ImageDirectionY;

  newScalarImage->SetImageOrigin( imageOrigin );
  newScalarImage->SetImageDirectionX( this->m_ImageDirectionX );
  newScalarImage->SetImageDirectionY( this->m_ImageDirectionY );
  newScalarImage->SetImageSlicePosition( this->m_ImageSlicePosition );
  
  newScalarImage->CreatePixelData( this->m_PixelData->GetType() );

  // compensate for dimensions that connot be evenly divided by downsampling 
  // factor.
  const unsigned int dimsY = (this->m_Dims[1] / factorY) * factorY;
  const unsigned int dimsX = (this->m_Dims[0] / factorX) * factorX;
  const Types::DataItem factorXY = 1.0 / (factorX * factorY);

  int j = 0;
  for ( unsigned int y = 0; y < dimsY; y += factorY, ++j ) 
    {
    int i = 0;
    for ( unsigned int x = 0; x < dimsX; x += factorX, ++i ) 
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

TypedArray* 
ScalarImage::GetMedianFiltered( const byte range ) const
{
  if ( !this->m_PixelData ) return NULL;

  TypedArray *result = this->m_PixelData->NewTemplateArray( this->m_PixelData->GetDataSize() );

  Types::DataItem *sort = Memory::AllocateArray<Types::DataItem>( range * range * range );

  unsigned delta = (range-1) / 2;
  unsigned offset = 0;
  for ( unsigned y = 0; y < this->m_Dims[1]; ++y )
    for ( unsigned x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
      {
      const unsigned int xFrom = ( x > delta ) ? ( x - delta ) : 0;
      const unsigned int yFrom = ( y > delta ) ? ( y - delta ) : 0;
      const unsigned int xTo = std::min( x+delta+1, this->m_Dims[0] );
      const unsigned int yTo = std::min( y+delta+1, this->m_Dims[1] );
      
      unsigned int source = 0;
      for ( unsigned int yy = yFrom; yy < yTo; ++yy )
	for ( unsigned int xx = xFrom; xx < xTo; ++xx, ++source ) 
	  {
	  this->m_PixelData->Get( sort[source], xx + this->m_Dims[0] * yy );
	  }
      
#ifdef CMTK_DATA_FLOAT	
      qsort( sort, source, sizeof( *sort ), MathUtil::CompareFloat );
#else
      qsort( sort, source, sizeof( *sort ), MathUtil::CompareDouble );
#endif

      if ( source % 2 )
	result->Set( sort[source/2], offset );
      else
	result->Set( (Types::DataItem) 
		     ( 0.5 * (sort[source/2] + sort[source/2-1]) ), offset );
    }

  delete[] sort;

  return result;
}

TypedArray*
ScalarImage::GetGaussFiltered( const Types::Coordinate stdDev ) const
{
  const Types::Coordinate stdDevPixelX = stdDev / this->m_PixelSize[0];
  const Types::Coordinate stdDevPixelY = stdDev / this->m_PixelSize[1];

  const size_t stdDevDiscreteX = static_cast<size_t>( ceil( stdDevPixelX ) );
  const size_t stdDevDiscreteY = static_cast<size_t>( ceil( stdDevPixelY ) );

  const size_t filterLengthX = std::min<size_t>( this->m_Dims[0], 3 * stdDevDiscreteX + 1 );
  const size_t filterLengthY = std::min<size_t>( this->m_Dims[1], 3 * stdDevDiscreteY + 1 );

  std::vector<Types::DataItem> filterX( filterLengthX );
  for ( size_t x=0; x < filterLengthX; ++x ) 
    {
    filterX[x] = 1.0/(sqrt(2*M_PI) * stdDevPixelX) * exp( -MathUtil::Square( 1.0 * x / stdDevPixelX ) / 2 );
    }
  
  std::vector<Types::DataItem> filterY( filterLengthY );
  for ( size_t y=0; y < filterLengthY; ++y ) 
    {
    filterY[y] = 1.0/(sqrt(2*M_PI) * stdDevPixelY) * exp( -MathUtil::Square( 1.0 * y / stdDevPixelY ) / 2);
    }
  
  TypedArray *result = this->GetFilteredData( filterX, filterY );
  
  return result;
}

TypedArray* 
ScalarImage::GetFilteredData
( const std::vector<Types::DataItem>& filterX, const std::vector<Types::DataItem>& filterY ) const
{
  if ( !this->m_PixelData ) return NULL;

  const size_t filterXsize = filterX.size();
  const size_t filterYsize = filterY.size();

  TypedArray *result = this->m_PixelData->NewTemplateArray( this->m_PixelData->GetDataSize() );

  size_t maxDim = std::max( this->m_Dims[0], this->m_Dims[1] );
  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );

  for ( unsigned int y=0; y < this->m_Dims[1]; ++y ) 
    {
    // copy row data to buffer
    for ( unsigned int x=0; x < this->m_Dims[0]; ++x )
      if ( !this->m_PixelData->Get( pixelBufferFrom[x], x+y*this->m_Dims[0] ) ) 
	pixelBufferFrom[x] = 0;
    
    // convolve row with filter
    for ( unsigned int x=0; x < this->m_Dims[0]; ++x ) 
      {
      // this keeps track of outside FOV data
      Types::DataItem correctOverlap = 0;
      // central element first to initialize target value
      pixelBufferTo[x] = pixelBufferFrom[x] * filterX[0];
      // now convolve side elements
      for ( unsigned int t=1; t < filterXsize; ++t ) 
	{
	// is x-t still inside the image?
	if ( x >= t )
	  // yes: convolve
	  pixelBufferTo[x] += pixelBufferFrom[x-t] * filterX[t];
	else
	  // no: save missing contribution for later
	  correctOverlap += filterX[t];

	// same for x+t:
	if ( x+t < this->m_Dims[0] )
	  pixelBufferTo[x] += pixelBufferFrom[x+t] * filterX[t];
	else
	  correctOverlap += filterX[t];
      }
      // correct value scaling for all missing (outside) elements
      pixelBufferTo[x] /= (1-correctOverlap);
    }
    
    for ( unsigned int x=0; x < this->m_Dims[0]; ++x )
      result->Set( pixelBufferTo[x], x+y*this->m_Dims[0] );
  }

  for ( unsigned int x=0; x < this->m_Dims[0]; ++x ) 
    {
    for ( unsigned int y=0; y < this->m_Dims[1]; ++y )
      if ( !result->Get( pixelBufferFrom[y], x+this->m_Dims[0]*y ) ) 
	pixelBufferFrom[y] = 0;

    for ( unsigned int y=0; y < this->m_Dims[1]; ++y ) 
      {
      Types::DataItem correctOverlap = 0;
      pixelBufferTo[y] = pixelBufferFrom[y] * filterY[0];
      for ( unsigned int t=1; t < filterYsize; ++t ) 
	{
	if ( y >= t )
	  pixelBufferTo[y] += pixelBufferFrom[y-t] * filterY[t];
	else
	  correctOverlap += filterY[t];
	
	if ( y+t < this->m_Dims[1] )
	  pixelBufferTo[y] += pixelBufferFrom[y+t] * filterY[t];
	else
	  correctOverlap += filterY[t];
      }
      // correct value scaling for all missing (outside) elements
      pixelBufferTo[y] /= (1-correctOverlap);
    }
    
    // write back convolved data
    for ( unsigned int y=0; y < this->m_Dims[1]; ++y )
      result->Set( pixelBufferTo[y], x+this->m_Dims[0]*y );
  }

  return result;
}

TypedArray*
ScalarImage::GetSobel2DFiltered() const
{
  if ( !this->m_PixelData ) return NULL;

  TypedArray *result = this->m_PixelData->NewTemplateArray( this->m_PixelData->GetDataSize() );
  
  Types::DataItem fov[3][3];
  size_t offset = 0;

  for ( unsigned int y = 0; y < this->m_Dims[1]; ++y )
    for ( unsigned int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
      {
      Types::DataItem value = 0;
      if ( x && y && (x<this->m_Dims[0]-1) && (y<this->m_Dims[1]-1) ) 
	{
	for ( int dy=0; dy<3; ++dy )
	  for ( int dx=0; dx<3; ++dx )
	    this->m_PixelData->Get( fov[dx][dy], x-1+dx+ this->m_Dims[0] * ( y-1+dy ) );
	
	value = 
	  sqrt( MathUtil::Square( fov[0][0] - fov[2][0] + 2 * ( fov[0][1] - fov[2][1] ) + fov[0][2] - fov[2][2] ) +
		MathUtil::Square( fov[0][0] - fov[0][2] + 2 * ( fov[1][0] - fov[1][2] ) + fov[2][0] - fov[2][2] ) );
	}
      result->Set( value, offset );
      }
  
  return result;
}

TypedArray* 
ScalarImage::GetSobelFiltered( const bool horizontal, const bool absolute ) 
  const
{
  if ( !this->m_PixelData ) return NULL;

  TypedArray *result = ( absolute ) ?
    TypedArray::Create( GetUnsignedDataType( this->m_PixelData->GetType() ), this->m_PixelData->GetDataSize() ) :
    TypedArray::Create( GetSignedDataType( this->m_PixelData->GetType() ), this->m_PixelData->GetDataSize() );
  
  Types::DataItem fov[3][3];
  size_t offset = 0;

  if ( horizontal )
    for ( unsigned int y = 0; y < this->m_Dims[1]; ++y )
      for ( unsigned int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
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
    for ( unsigned int y = 0; y < this->m_Dims[1]; ++y )
      for ( unsigned int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
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

TypedArray*
ScalarImage::GetLaplace2DFiltered() const
{
  if ( this->m_PixelData.IsNull() ) return NULL;

  TypedArray *result = this->m_PixelData->NewTemplateArray( this->m_PixelData->GetDataSize() );
  
  Types::DataItem fov[3][3];
  size_t offset = 0;

  for ( unsigned int y = 0; y < this->m_Dims[1]; ++y )
    for ( unsigned int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
      {
      Types::DataItem value = 0;
      if ( x && y && (x<this->m_Dims[0]-1) && (y<this->m_Dims[1]-1) ) 
	{
	for ( int dy=0; dy<3; ++dy )
	  for ( int dx=0; dx<3; ++dx )
	    this->m_PixelData->Get( fov[dx][dy], x-1+dx+ this->m_Dims[0] * ( y-1+dy ) );
	
	value = fov[0][1] + fov[2][1] + fov[1][0] + fov[1][2] + fov[0][0] + fov[0][2] + fov[2][0] + fov[2][2] - 8 * fov[1][1];
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
( const ScalarImage& other )
{
  TypedArray *data0 = this->GetPixelData().GetPtr();
  const TypedArray *data1 = other.GetPixelData();

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
    for ( unsigned int y = 0; y < this->m_Dims[1]/2; ++y ) 
      {
      this->m_PixelData->BlockSwap( y * this->m_Dims[0], (this->m_Dims[1]-y-1) * this->m_Dims[0], this->m_Dims[0] );
      }
    this->m_ImageOrigin = this->m_ImageOrigin + ((this->m_Dims[1]-1) * this->m_PixelSize[1] / this->m_ImageDirectionY.EuclidNorm()) * this->m_ImageDirectionY;
    this->m_ImageDirectionY *= (-1.0);
    }
  
  if ( horizontal ) 
    {
    for ( unsigned int y = 0; y < this->m_Dims[1]; ++y ) 
      {
      this->m_PixelData->BlockReverse( y * this->m_Dims[0], this->m_Dims[0] );
      }
    this->m_ImageOrigin = this->m_ImageOrigin + ((this->m_Dims[1]-1) * this->m_PixelSize[0] / this->m_ImageDirectionX.EuclidNorm()) * this->m_ImageDirectionX;
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
  unsigned int newDimsY = static_cast<unsigned int>( (this->m_Dims[1]-1) * this->m_PixelSize[1]/this->m_PixelSize[0] ) + 1;
  
  TypedArray *scaled = this->m_PixelData->NewTemplateArray( this->m_Dims[0] * newDimsY );
  
  if ( interpolate ) 
    {
    // with interpolation
    std::vector<Types::DataItem> row0( this->m_Dims[0] ), row1( this->m_Dims[0] );
    this->m_PixelData->GetSubArray( &(row0[0]), 0, this->m_Dims[0] );
    this->m_PixelData->GetSubArray( &(row1[0]), this->m_Dims[0], this->m_Dims[0] );
    
    Types::Coordinate scanLine = 0;
    size_t ySource = 0;
    size_t offset = 0;
    for ( unsigned int y = 0; y < newDimsY; ++y ) 
      {
      Types::Coordinate factor = scanLine / this->m_PixelSize[1];
      for ( unsigned int x = 0; x < this->m_Dims[0]; ++x, ++offset ) 
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
    size_t ySource = 0;
    for ( unsigned int y = 0; y < newDimsY; ++y ) 
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
  this->SetPixelData( TypedArray::SmartPtr( scaled ) );
}

void
ScalarImage::AdjustAspectX( const bool interpolate )
{
  unsigned int newDimsX = static_cast<unsigned int>( (this->m_Dims[0]-1) * this->m_PixelSize[0]/this->m_PixelSize[1] ) + 1;
  
  TypedArray *scaled = this->m_PixelData->NewTemplateArray( newDimsX * this->m_Dims[1] );
  
  if ( interpolate ) 
    {
    std::vector<Types::Coordinate> factor( newDimsX );
    std::vector<unsigned int> fromPixel( newDimsX );
    
    Types::Coordinate scanLine = 0;
    size_t xSource = 0;
    for ( unsigned int x = 0; x < newDimsX; ++x ) 
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
    for ( unsigned int y = 0; y < this->m_Dims[1]; ++y ) 
      {
      this->m_PixelData->GetSubArray( &(rowFrom[0]), y * this->m_Dims[0], this->m_Dims[0] );
      for ( unsigned int x = 0; x < newDimsX; ++x, ++offset ) 
	{
	scaled->Set( (1.0 - factor[x] ) * rowFrom[fromPixel[x]] + factor[x] * rowFrom[fromPixel[x]+1], offset );
	}
      }
    } 
  else
    {
    Types::Coordinate scanLine = this->m_PixelSize[0] / 2; // correct offset for NN
    size_t xSource = 0;
    std::vector<unsigned int> fromPixel( newDimsX );
    for ( unsigned int x = 0; x < newDimsX; ++x ) 
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
    
    for ( unsigned int y = 0; y < this->m_Dims[1]; ++y ) 
      {
      for ( unsigned int x = 0; x < newDimsX; ++x ) 
	{
	memcpy( scaledPtr, pixelPtr+fromPixel[x], scaled->GetItemSize() );
	scaledPtr += scaled->GetItemSize();
	}
      pixelPtr += scaled->GetItemSize() * this->m_Dims[0];
      }
    }
  
  this->m_PixelSize[0] = this->m_PixelSize[1];
  this->m_Dims[0] = newDimsX;
  this->SetPixelData( TypedArray::SmartPtr( scaled ) );
}

void
ScalarImage::ProjectPixel
( const Vector3D& v, unsigned int& i, unsigned int& j ) const
{
  Vector3D p(v);
  p.XYZ[0] -= this->m_ImageOrigin[0];
  p.XYZ[1] -= this->m_ImageOrigin[1];
  p.XYZ[2] -= this->m_ImageOrigin[2];

  i = MathUtil::Round( ( p * this->m_ImageDirectionX ) / ( this->m_ImageDirectionX.Square() * this->m_PixelSize[0] ) );
  j = MathUtil::Round( ( p * this->m_ImageDirectionY ) / ( this->m_ImageDirectionY.Square() * this->m_PixelSize[1] ) );
}

void 
ScalarImage::Print() const
{
  StdErr.printf( "ScalarImage at %p\n", this );

  StdErr.indent();
  StdErr.printf( "Dimensions: [%d,%d]\n", this->m_Dims[0], this->m_Dims[1] );
  StdErr.printf( "Pixel size: [%f,%f]\n", this->m_PixelSize[0], this->m_PixelSize[1] );
  StdErr.printf( "Origin: [%f,%f,%f]\n", this->m_ImageOrigin.XYZ[0], this->m_ImageOrigin.XYZ[1], this->m_ImageOrigin.XYZ[2] );

  StdErr.printf( "DirectionX: [%f,%f,%f]\n", this->m_ImageDirectionX.XYZ[0], this->m_ImageDirectionX.XYZ[1], this->m_ImageDirectionX.XYZ[2] );
  StdErr.printf( "DirectionY: [%f,%f,%f]\n", this->m_ImageDirectionY.XYZ[0], this->m_ImageDirectionY.XYZ[1], this->m_ImageDirectionY.XYZ[2] );
  StdErr.unindent();
}

} // namespace cmtk
