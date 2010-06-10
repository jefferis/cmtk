/*
//
//  Copyright 1997-2009 Torsten Rohlfing
//
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

#include <cmtkDataGridMorphologicalOperators.h>

#include <cmtkTemplateArray.h>
#include <cmtkException.h>

namespace
cmtk
{

DataGridMorphologicalOperators
::DataGridMorphologicalOperators( const DataGrid::SmartConstPtr& dataGrid ) 
  : m_DataGrid( dataGrid ) 
{
  if ( !this->m_DataGrid->GetData() )
    throw Exception( "ERROR: DataGrid object given to DataGridMorphologicalOperators constructor does not have a data array" );
}

bool
DataGridMorphologicalOperators::EliminatePaddingVoting( const int iterations )
{
  bool changed = false;
  for ( int it = 0; (it < iterations) && changed; ++it )
    {
    const TypedArray* data = this->m_DataGrid->GetData();
    TypedArray* result = data->Clone();
    
    if ( !data->GetPaddingFlag() ) return changed;
    const Types::DataItem dataPadding = data->GetPaddingValue();
    
    short valueMap[256];  
    size_t offset = 0;
    for ( int z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
      {
      const int dzFrom = z ? z-1 : z, dzTo = (z<this->m_DataGrid->m_Dims[2]-1) ? z+1 : z;
      for ( int y = 0; y < this->m_DataGrid->m_Dims[1]; ++y ) 
	{
	const int dyFrom = y ? y-1 : y, dyTo = (y<this->m_DataGrid->m_Dims[1]-1) ? y+1 : y;
	for ( int x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	  {
	  // is there Padding here?
	  if ( data->PaddingDataAt( offset ) )
	    {
	    // determine frequencies of values among neighbors
	    memset( valueMap, 0, sizeof( valueMap ) );
	    const int dxFrom = x ? x-1 : x, dxTo = (x<this->m_DataGrid->m_Dims[0]-1) ? x+1 : x;
	    for ( int dz = dzFrom; (dz <= dzTo); ++dz )
	      for ( int dy = dyFrom; (dy <= dyTo); ++dy )
		for ( int dx = dxFrom; (dx <= dxTo); ++dx )
		  {
		  Types::DataItem value;
		  if ( this->m_DataGrid->GetDataAt( value, dx, dy, dz ) )
		    {
		    ++valueMap[static_cast<byte>( value )];
		    }
		  }
	    
	    // find most frequent value among neighbors
	    char maxCount = 0;
	    Types::DataItem maxValue = 0;
	    
	    for ( int it = 0; it < 256; ++it )
	      {
	      if ( valueMap[it] > maxCount )
		{
		maxCount = valueMap[it];
		maxValue = it;
		}
	      else
		if ( valueMap[it] == maxCount )
		  {
		  maxValue = dataPadding;
		  }
	      }
	    
	    if ( maxValue != dataPadding )
	      {
	      changed = true;
	      result->Set( maxValue, offset );
	      }
	    } // if PaddingDataAt()
	  } // for x
	} // for y
      } // for z
    }
  
  return changed;
}

TypedArray::SmartPtr
DataGridMorphologicalOperators::GetEroded( const int iterations ) const
{
  TypedArray::SmartPtr dataArray = this->m_DataGrid->GetData();
  if ( dataArray->GetType() != TYPE_BYTE ) 
    {
    throw Exception( "ERROR: convert data to byte before calling DataGridMorphologicalOperators::GetDataErode()" );
    }

  const byte* data = static_cast<const byte*>( dataArray->GetDataPtr() );
  byte *tmp = Memory::AllocateArray<byte>(  dataArray->GetDataSize()  );

  ByteArray::SmartPtr erodedArray = ByteArray::Create( dataArray->GetDataSize() );
  byte* eroded = erodedArray->GetDataPtrConcrete();
  
  memcpy( eroded, data, erodedArray->GetDataSizeBytes() );
  
  for ( int i = 0; i < iterations; ++i ) 
    {
    size_t offset = 0;
    for ( int z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
      {
      int dzFrom = z ? -1 : 0, dzTo = (z<this->m_DataGrid->m_Dims[2]-1) ? 1 : 0;
      for ( int y = 0; y < this->m_DataGrid->m_Dims[1]; ++y ) 
	{
	int dyFrom = y ? -1 : 0, dyTo = (y<this->m_DataGrid->m_Dims[1]-1) ? 1 : 0;
	for ( int x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	  {
	  int dxFrom = x ? -1 : 0, dxTo = (x<this->m_DataGrid->m_Dims[0]-1) ? 1 : 0;
	  if ( eroded[offset] ) 
	    {
	    bool erodePixel = false;
	    for ( int dz = dzFrom; (dz <= dzTo) && !erodePixel; ++dz )
	      for ( int dy = dyFrom; (dy <= dyTo) && !erodePixel; ++dy )
		for ( int dx = dxFrom; (dx <= dxTo) && !erodePixel; ++dx )
		  if ( dx || dy || dz )
		    if ( ! eroded[offset+this->m_DataGrid->GetOffsetFromIndex( dx, dy, dz )] )
		      erodePixel = true;
	    if ( erodePixel )
	      tmp[offset] = 0;
	    else
	      tmp[offset] = eroded[offset];
	    } 
	  else
	    { // if data
	    tmp[offset] = 0;
	    } // if data
	  } // for x
	} // for y
      } // for z
    
    memcpy( eroded, tmp, erodedArray->GetDataSizeBytes() );
    } // for i
  
  delete[] tmp;
  
  return erodedArray;
}

TypedArray::SmartPtr
DataGridMorphologicalOperators::GetDilated( const int iterations ) const
{
  TypedArray::SmartPtr dataArray = this->m_DataGrid->GetData();
  if ( dataArray->GetType() != TYPE_BYTE ) 
    {
    throw Exception( "ERROR: convert data to byte before calling DataGridMorphologicalOperators::GetDataDilate()" );
    }
  
  const byte* data = static_cast<const byte*>( dataArray->GetDataPtr() );
  byte *tmp = Memory::AllocateArray<byte>(  dataArray->GetDataSize()  );
  
  ByteArray::SmartPtr dilatedArray = ByteArray::Create( dataArray->GetDataSize() );
  byte* dilated = dilatedArray->GetDataPtrConcrete();
  
  memcpy( dilated, data, dilatedArray->GetDataSizeBytes() );
  
  for ( int i = 0; i < iterations; ++i ) 
    {
    size_t offset = 0;
    for ( int z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
      {
      int dzFrom = z ? -1 : 0, dzTo = (z<this->m_DataGrid->m_Dims[2]-1) ? 1 : 0;
      for ( int y = 0; y < this->m_DataGrid->m_Dims[1]; ++y ) 
	{
	int dyFrom = y ? -1 : 0, dyTo = (y<this->m_DataGrid->m_Dims[1]-1) ? 1 : 0;
	for ( int x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	  {
	  int dxFrom = x ? -1 : 0, dxTo = (x<this->m_DataGrid->m_Dims[0]-1) ? 1 : 0;
	  if ( ! dilated[offset] ) 
	    {
	    byte dilatePixel = 0;
	    for ( int dz = dzFrom; (dz <= dzTo) && !dilatePixel; ++dz )
	      for ( int dy = dyFrom; (dy <= dyTo) && !dilatePixel; ++dy )
		for ( int dx = dxFrom; (dx <= dxTo) && !dilatePixel; ++dx )
		  if ( dx || dy || dz )
		    dilatePixel = dilated[offset+this->m_DataGrid->GetOffsetFromIndex(dx,dy,dz)];
	    if ( dilatePixel )
	      tmp[offset] = dilatePixel;
	    else
	      tmp[offset] = 0;
	    } 
	  else
	    { // !data
	    tmp[offset] = dilated[offset];
	    } // else !data
	  } // for x
	} // for y
      } // for z
    
    memcpy( dilated, tmp, dilatedArray->GetDataSizeBytes() );
    } // for i
  
  delete[] tmp;
  
  return dilatedArray;
}

TypedArray::SmartPtr 
DataGridMorphologicalOperators::GetBoundaryMap( const bool multiValued ) const
{
  TypedArray::SmartPtr dataArray = this->m_DataGrid->GetData();
  ShortArray::SmartPtr boundaryArray = ShortArray::Create( dataArray->GetDataSize() );
  short* boundary = boundaryArray->GetDataPtrConcrete();

#pragma omp parallel for
  for ( int z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
    {
    size_t offset = z * this->m_DataGrid->m_Dims[0] * this->m_DataGrid->m_Dims[1];

    Types::DataItem value, neighbor;
    const int dzFrom = z ? -1 : 0, dzTo = (z<this->m_DataGrid->m_Dims[2]-1) ? 1 : 0;
    for ( int y = 0; y < this->m_DataGrid->m_Dims[1]; ++y ) 
      {
      const int dyFrom = y ? -1 : 0, dyTo = (y<this->m_DataGrid->m_Dims[1]-1) ? 1 : 0;
      for ( int x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	{
	const int dxFrom = x ? -1 : 0, dxTo = (x<this->m_DataGrid->m_Dims[0]-1) ? 1 : 0;
	bool bp = false;
	if ( dataArray->Get( value, offset ) ) 
	  {
	  for ( int dz = dzFrom; (dz <= dzTo) && !bp; ++dz )
	    for ( int dy = dyFrom; (dy <= dyTo) && !bp; ++dy )
	      for ( int dx = dxFrom; (dx <= dxTo) && !bp; ++dx )
		if ( dx || dy || dz )
		  if ( dataArray->Get( neighbor, offset+this->m_DataGrid->GetOffsetFromIndex(dx,dy,dz) ) )
		    bp = (value != neighbor);
	  }
	else
	  {
	  value = 0;
	  }

	if ( bp )
	  {
	  if ( multiValued )
	    boundary[offset] = static_cast<short>( value );
	  else
	    boundary[offset] = 1;
	  }
	else
	  boundary[offset] = 0;
	
	} // for x
      } // for y
    } // for z
  
  return boundaryArray;
}

} // namespace cmtk
