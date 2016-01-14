/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2011, 2013-2014 SRI International
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

#include <algorithm>

template<class TFilter>
cmtk::TypedArray::SmartPtr
cmtk::DataGridFilter::ApplyRegionFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  const TypedArray* data = this->m_DataGrid->GetData();
  if ( !data )
    return TypedArray::SmartPtr( NULL );

  TypedArray::SmartPtr result = TypedArray::Create( data->GetType(), data->GetDataSize() );
  
  const Types::GridIndexType widthX = 1 + 2*radiusX;
  const Types::GridIndexType widthY = 1 + 2*radiusY;
  const Types::GridIndexType widthZ = 1 + 2*radiusZ;

  const Types::GridIndexType pixelsPerPlane = this->m_DataGrid->m_Dims[0] * this->m_DataGrid->m_Dims[1];
#pragma omp parallel for
  for ( Types::GridIndexType z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
    {
    Types::GridIndexType offset = z * pixelsPerPlane;
    std::vector<Types::DataItem> regionData(  widthX*widthY*widthZ  );
  
    Types::GridIndexType zFrom = ( z > radiusZ ) ? ( z - radiusZ ) : 0;
    Types::GridIndexType zTo = std::min( z+radiusZ+1, this->m_DataGrid->m_Dims[2] );
    
    for ( Types::GridIndexType y = 0; y < this->m_DataGrid->m_Dims[1]; ++y ) 
      {      
      Types::GridIndexType yFrom = ( y > radiusY ) ? ( y - radiusY ) : 0;
      Types::GridIndexType yTo = std::min( y+radiusY+1, this->m_DataGrid->m_Dims[1] );
      
      for ( Types::GridIndexType x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	{
	Types::GridIndexType xFrom = ( x > radiusX ) ? ( x - radiusX ) : 0;
	Types::GridIndexType xTo = std::min( x+radiusX+1, this->m_DataGrid->m_Dims[0] );
	
	regionData.clear();

	Types::GridIndexType ofsZ = yFrom + this->m_DataGrid->m_Dims[1] * zFrom;
	for ( Types::GridIndexType zz = zFrom; zz < zTo; ++zz, ofsZ += this->m_DataGrid->m_Dims[1] ) 
	  {
	  Types::GridIndexType ofsYZ = this->m_DataGrid->m_Dims[0] * ofsZ ;
	  for ( Types::GridIndexType yy = yFrom; yy < yTo; ++yy, ofsYZ += this->m_DataGrid->m_Dims[0] ) 
	    {
	    Types::GridIndexType toYZ = ofsYZ + xTo;
	    for ( Types::GridIndexType xx = xFrom + ofsYZ; xx < toYZ; ++xx ) 
	      {
	      Types::DataItem value = 0;
	      if ( data->Get( value, xx ) )
		{
		regionData.push_back( value );
		}
	      }
	    }
	  }

	result->Set( TFilter::Reduce( regionData ), offset );
	}
      }
    }
  return result;
}

