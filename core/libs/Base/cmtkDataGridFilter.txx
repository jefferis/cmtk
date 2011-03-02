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

#include <System/cmtkProgress.h>

#include <algorithm>

template<class TFilter>
cmtk::TypedArray::SmartPtr
cmtk::DataGridFilter::ApplyRegionFilter( const int radiusX, const int radiusY, const int radiusZ ) const
{
  const TypedArray* data = this->m_DataGrid->GetData();
  TypedArray::SmartPtr result = TypedArray::Create( data->GetType(), data->GetDataSize() );
  
  const int widthX = 1 + 2*radiusX;
  const int widthY = 1 + 2*radiusY;
  const int widthZ = 1 + 2*radiusZ;

  std::vector<Types::DataItem> regionData(  widthX*widthY*widthZ  );
  
  int offset = 0;
  Progress::Begin( 0, this->m_DataGrid->m_Dims[2], 1 );

  Progress::ResultEnum status = Progress::OK;
  for ( int z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
    {
    status = Progress::SetProgress( z );
    if ( status != Progress::OK ) break;
    
    int zFrom = ( z > radiusZ ) ? ( z - radiusZ ) : 0;
    int zTo = std::min( z+radiusZ+1, this->m_DataGrid->m_Dims[2] );
    
    for ( int y = 0; y < this->m_DataGrid->m_Dims[1]; ++y ) 
      {      
      int yFrom = ( y > radiusY ) ? ( y - radiusY ) : 0;
      int yTo = std::min( y+radiusY+1, this->m_DataGrid->m_Dims[1] );
      
      for ( int x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	{
	int xFrom = ( x > radiusX ) ? ( x - radiusX ) : 0;
	int xTo = std::min( x+radiusX+1, this->m_DataGrid->m_Dims[0] );
	
	regionData.clear();

	int ofsZ = yFrom + this->m_DataGrid->m_Dims[1] * zFrom;
	for ( int zz = zFrom; zz < zTo; ++zz, ofsZ += this->m_DataGrid->m_Dims[1] ) 
	  {
	  int ofsYZ = this->m_DataGrid->m_Dims[0] * ofsZ ;
	  for ( int yy = yFrom; yy < yTo; ++yy, ofsYZ += this->m_DataGrid->m_Dims[0] ) 
	    {
	    int toYZ = ofsYZ + xTo;
	    for ( int xx = xFrom + ofsYZ; xx < toYZ; ++xx ) 
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
  Progress::Done();
  
  if ( status != Progress::OK ) 
    {
    result = TypedArray::SmartPtr( NULL );
    }
  
  return result;
}

