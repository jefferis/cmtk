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

#include <cmtkImage.h>

namespace
cmtk
{

/** \addtogroup Pipeline */
//@{

Image::Image () :
  Data( NULL )
{
  DataType = TYPE_NONE;
}

TypedArray::SmartPtr
Image::GetData()
{
  if ( ! Data ) 
    {
    if ( DataType == TYPE_NONE )
      return TypedArray::SmartPtr( NULL );
    else
      {
      Data = TypedArray::SmartPtr( TypedArray::Create( DataType, Dims[0] * Dims[1] ) );
      this->UpdateModifiedTime();
      }
    } 
  else
    {
    if ( ( Data->GetType() != DataType ) || ( Data->GetDataSize() != (Dims[0] * Dims[1]) ) ) 
      {
      Data = TypedArray::SmartPtr( NULL );
      this->UpdateModifiedTime();
      return this->GetData();
      } 
    }
  return Data;
}

void
Image::SetData( TypedArray::SmartPtr& data )
{
  Data = data;
  if ( Data ) 
    DataType = Data->GetType();
  this->UpdateModifiedTime();
}

void
Image::SetFromScalarImage
( ScalarImage *const scalarImage, const bool copyPixelData )
{
  if ( scalarImage ) 
    {
    this->SetDims( scalarImage->GetDims()[0], scalarImage->GetDims()[1] );
    TypedArray::SmartPtr pixelData = scalarImage->GetPixelData();
    if ( copyPixelData )
      pixelData = TypedArray::SmartPtr( pixelData->Clone() );
    this->SetData( pixelData );
    this->SetSpacing( scalarImage->GetPixelSize() );
    this->SetOrigin( scalarImage->GetImageOrigin().begin() );
    this->SetDirectionX( scalarImage->GetImageDirectionX().begin() );
    this->SetDirectionY( scalarImage->GetImageDirectionY().begin() );
    this->UpdateModifiedTime();
    }
}

void
Image::SetFromScalarImage
( const ScalarImage* scalarImage )
{
  if ( scalarImage )
    {
    this->SetDims( scalarImage->GetDims()[0], scalarImage->GetDims()[1] );
    TypedArray::SmartPtr pixelData = scalarImage->GetPixelData();
    if ( pixelData )
      pixelData = TypedArray::SmartPtr( pixelData->Clone() );
    this->SetData( pixelData );
    this->SetSpacing( scalarImage->GetPixelSize() );
    this->SetOrigin( scalarImage->GetImageOrigin().begin() );
    this->SetDirectionX( scalarImage->GetImageDirectionX().begin() );
    this->SetDirectionY( scalarImage->GetImageDirectionY().begin() );
    this->UpdateModifiedTime();
    }
}

ScalarImage* 
Image::GetScalarImage() const
{
  ScalarImage* scalarImage = new ScalarImage( Dims[0], Dims[1] );

  scalarImage->SetPixelData( Data );
  scalarImage->SetPixelSize( Spacing );
  scalarImage->SetImageOrigin( FixedVector<3,Types::Coordinate>( Origin ) );
  scalarImage->SetImageDirectionX( FixedVector<3,Types::Coordinate>( DirectionX ) );
  scalarImage->SetImageDirectionY( FixedVector<3,Types::Coordinate>( DirectionY ) );

  return scalarImage;
}

double
Image::GetDataAt( const int x, const int y, const double def )
{
  const TypedArray* data = this->GetData();

  Types::DataItem result;
  if ( data->Get( result, x+Dims[0]*y ) ) 
    {
    return result;
    } 
  else
    {
    return def;
    }
}

double
Image::GetDataAt( const int index, const double def )
{
  const TypedArray *data = this->GetData();

  Types::DataItem result;
  if ( data->Get( result, index ) ) 
    {
    return result;
    } 
  else
    {
    return def;
    }
}

void
Image::SetDataAt( const int x, const int y, const double value )
{
  this->GetData()->Set( value, x+Dims[0]*y );
}

void
Image::SetDataAt( const int index, const double value )
{
  this->GetData()->Set( value, index );
}

double
Image::GetDataAt( const double x, const double y, const double def )
{
  const TypedArray *data = this->GetData();

  const unsigned int idxX = static_cast<int>( x / Spacing[0] );
  const unsigned int idxY = static_cast<int>( y / Spacing[1] );

  if ( (idxX > Dims[0]-2) || (idxY > Dims[1]-2) )
    return def;

  int offset = idxX + Dims[0] * idxY;
  
  Types::DataItem result[4];
  if ( ! data->Get( result[0], offset ) ) 
    {
    return def;
    }
  if ( ! data->Get( result[1], offset + 1 ) ) 
    {
    return def;
    }
  if ( ! data->Get( result[2], offset + Dims[0] ) ) 
    {
    return def;
    }
  if ( ! data->Get( result[3], offset + Dims[0] + 1 ) ) 
    {
    return def;
    }
  
  const double relX = ( x - idxX * Spacing[0] ) / Spacing[0];
  const double relY = ( y - idxY * Spacing[1] ) / Spacing[1];

  return (1-relY) * ( (1-relX) * result[0] + (relX) * result[1] ) + (relY) * ( (1-relX) * result[2] + (relX) * result[3] );
}

} // namespace cmtk
