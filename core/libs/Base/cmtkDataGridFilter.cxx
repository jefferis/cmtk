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

#include <cmtkDataGridFilter.h>

#include <cmtkException.h>
#include <cmtkMemory.h>
#include <cmtkProgress.h>
#include <cmtkThreadPool.h>

#include <vector>

namespace
cmtk
{

DataGridFilter
::DataGridFilter( DataGrid::SmartPtr dataGrid )
  : m_DataGrid( dataGrid )
{
  if ( !this->m_DataGrid->GetData() )
    throw Exception( "ERROR: DataGrid object given to DataGridFilter constructor does not have a data array" );
}

TypedArray::SmartPtr
DataGridFilter::GetDataMedianFiltered( const int radiusX, const int radiusY, const int radiusZ ) const
{
  const TypedArray* data = this->m_DataGrid->GetData();
  TypedArray::SmartPtr result = TypedArray::Create( data->GetType(), data->GetDataSize() );
  
  const int widthX = 1 + 2*radiusX;
  const int widthY = 1 + 2*radiusY;
  const int widthZ = 1 + 2*radiusZ;
  Types::DataItem *sort = Memory::AllocateArray<Types::DataItem>(  widthX*widthY*widthZ  );
  
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
	
	int source = 0;
	int ofsZ = yFrom + this->m_DataGrid->m_Dims[1] * zFrom;
	for ( int zz = zFrom; zz < zTo; ++zz, ofsZ += this->m_DataGrid->m_Dims[1] ) 
	  {
	  int ofsYZ = this->m_DataGrid->m_Dims[0] * ofsZ ;
	  for ( int yy = yFrom; yy < yTo; ++yy, ofsYZ += this->m_DataGrid->m_Dims[0] ) 
	    {
	    int toYZ = ofsYZ + xTo;
	    for ( int xx = xFrom + ofsYZ; xx < toYZ; ++xx, ++source ) 
	      {
	      data->Get( sort[source], xx );
	      }
	    }
	  }
	
#ifdef CMTK_DATA_FLOAT	
	qsort( sort, source, sizeof( *sort ), MathUtil::CompareFloat );
#else
	qsort( sort, source, sizeof( *sort ), MathUtil::CompareDouble );
#endif
	
	if ( source % 2 )
	  result->Set( sort[source/2], offset );
	else
	  result->Set( (Types::DataItem) ( 0.5 * (sort[source/2] + sort[source/2-1]) ), offset );
	}
      }
    }
  Progress::Done();
  
  delete[] sort;
  
  if ( status != Progress::OK ) 
    {
    result = TypedArray::SmartPtr( NULL );
    }
  
  return result;
}

TypedArray::SmartPtr
DataGridFilter::GetDataSobelFiltered() const
{
  const TypedArray* data = this->m_DataGrid->GetData();
  TypedArray::SmartPtr result = TypedArray::Create( data->GetType(), data->GetDataSize() );

  Types::DataItem value = 0;
  Types::DataItem fov[3][3][3];

  Progress::Begin( 0, this->m_DataGrid->m_Dims[2], 1 );

  size_t offset = 0;
  for ( int z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
    {
    Progress::SetProgress( z );
    for ( int y = 0; y < this->m_DataGrid->m_Dims[1]; ++y )
      for ( int x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	{
	if ( x && y && z && ( x<this->m_DataGrid->m_Dims[0]-1 ) && ( y<this->m_DataGrid->m_Dims[1]-1 ) && ( z<this->m_DataGrid->m_Dims[2]-1 ) ) 
	  {
	  for ( int dz=-1; dz<2; ++dz )
	    for ( int dy=-1; dy<2; ++dy )
	      for ( int dx=-1; dx<2; ++dx )
		if ( ! data->Get( fov[1+dx][1+dy][1+dz], offset+this->m_DataGrid->GetOffsetFromIndex( dx, dy, dz ) ) )
		  fov[1+dx][1+dy][1+dz] = 0;
	  
	  value = (Types::DataItem)
	    ( fabs( fov[0][0][1] - fov[2][0][1] + 
		    2 * ( fov[0][1][1] - fov[2][1][1] ) +
		    fov[0][2][1] - fov[2][2][1] ) +
	      fabs( fov[0][0][1] - fov[0][2][1] + 
		    2 * ( fov[1][0][1] - fov[1][2][1] ) +
		    fov[2][0][1] - fov[2][2][1] )+
	      fabs( fov[0][1][0] - fov[2][1][0] + 
		    2 * ( fov[0][1][1] - fov[2][1][1] ) +
		    fov[0][1][2] - fov[2][1][2] ) +
	      fabs( fov[0][1][0] - fov[0][1][2] + 
		    2 * ( fov[1][1][0] - fov[1][1][2] ) +
		    fov[2][1][0] - fov[2][1][2] ) +
	      fabs( fov[1][0][0] - fov[1][2][0] + 
		    2 * ( fov[1][0][1] - fov[1][2][1] ) +
		    fov[1][0][2] - fov[1][2][2] ) +
	      fabs( fov[1][0][0] - fov[1][0][2] + 
		    2 * ( fov[1][1][0] - fov[1][1][2] ) +
		    fov[1][2][0] - fov[1][2][2] ) ) / 6;
	  
	  result->Set( value, offset );
	  } 
	else
	  {
	  result->Set( 0, offset );
	  }
	}
    }
  
  Progress::Done();

  return result;
}

TypedArray::SmartPtr
DataGridFilter::GetDataKernelFiltered
( const std::vector<Types::DataItem>& filterX, 
  const std::vector<Types::DataItem>& filterY,
  const std::vector<Types::DataItem>& filterZ ) const
{
  TypedArray::SmartPtr result = TypedArray::Create( this->m_DataGrid->GetData()->GetType(), this->m_DataGrid->GetNumberOfPixels() );
  
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const size_t numberOfTasks = 4 * threadPool.GetNumberOfThreads() - 3;

  std::vector<FilterThreadParameters> params( numberOfTasks );
  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    params[task].thisObject = this;
    params[task].m_Filter = &filterX;
    params[task].m_Result = result;
    }
  threadPool.Run( GetFilteredDataThreadX, params );

  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    params[task].m_Filter = &filterY;
    }
  threadPool.Run( GetFilteredDataThreadY, params );

  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    params[task].m_Filter = &filterZ;
    }
  threadPool.Run( GetFilteredDataThreadZ, params );
  
  return result;
}

void
DataGridFilter
::GetFilteredDataThreadX( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  FilterThreadParameters* params = static_cast<FilterThreadParameters*>( args );
  const Self* ThisConst = params->thisObject;
  
  const DataGrid::IndexType& dims = ThisConst->m_DataGrid->m_Dims;
  unsigned int maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const size_t filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray::SmartPtr& result = params->m_Result;

  for ( int z=taskIdx; z < dims[2]; z += taskCnt ) 
    {
    for ( int y=0; y < dims[1]; ++y ) 
      {
      // copy row data to buffer
      size_t ofs = ThisConst->m_DataGrid->GetOffsetFromIndex( 0, y, z );
      for ( int x=0; x < dims[0]; ++x, ++ofs )
	if ( !ThisConst->m_DataGrid->GetDataAt( pixelBufferFrom[x], ofs ) )
	  pixelBufferFrom[x] = 0;
      
      // convolve row with filter
      for ( int x=0; x < dims[0]; ++x ) 
	{
	// this keeps track of outside FOV data
	Types::DataItem correctOverlap = filter[0];
	// central element first to initialize target value
	pixelBufferTo[x] = pixelBufferFrom[x] * filter[0];
	// now convolve side elements
	for ( int t=1; t < (int)filterSize; ++t ) 
	  {
	  // is x-t still inside the image?
	  if ( x >= t )
	    {
	    // yes: convolve
	    pixelBufferTo[x] += pixelBufferFrom[x-t] * filter[t];
	    correctOverlap += filter[t];
	    }
	  
	  // same for x+t:
	  if ( x+t < dims[0] )
	    {
	    pixelBufferTo[x] += pixelBufferFrom[x+t] * filter[t];
	    correctOverlap += filter[t];
	    }
	  }
	// correct value scaling for all missing (outside) elements
	pixelBufferTo[x] /= correctOverlap;
	}
      
      ofs = ThisConst->m_DataGrid->GetOffsetFromIndex( 0, y, z );
      for ( int x=0; x < dims[0]; ++x, ++ofs )
	result->Set( pixelBufferTo[x], ofs );
      }
    }
}

void
DataGridFilter
::GetFilteredDataThreadY( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  FilterThreadParameters* params = static_cast<FilterThreadParameters*>( args );
  const Self* ThisConst = params->thisObject;
  
  const DataGrid::IndexType& dims = ThisConst->m_DataGrid->m_Dims;
  unsigned int maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const size_t filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray::SmartPtr& result = params->m_Result;

  for ( int z=taskIdx; z < dims[2]; z+=taskCnt ) 
    {
    for ( int x=0; x < dims[0]; ++x ) 
      {
      for ( int y=0; y < dims[1]; ++y )
	if ( !result->Get( pixelBufferFrom[y], ThisConst->m_DataGrid->GetOffsetFromIndex( x, y, z ) ) ) 
	  pixelBufferFrom[y] = 0;
      
      for ( int y=0; y < dims[1]; ++y ) 
	{
	Types::DataItem correctOverlap = filter[0];
	pixelBufferTo[y] = pixelBufferFrom[y] * filter[0];
	for ( int t=1; t < (int)filterSize; ++t ) 
	  {
	  if ( y >= t ) 
	    {
	    pixelBufferTo[y] += pixelBufferFrom[y-t] * filter[t];
	    correctOverlap += filter[t];
	    }
	  
	  if ( y+t < dims[1] )
	    {
	    pixelBufferTo[y] += pixelBufferFrom[y+t] * filter[t];
	    correctOverlap += filter[t];
	    }
	  }
	// correct value scaling for all missing (outside) elements
	pixelBufferTo[y] /= correctOverlap;
	}
      
      // write back convolved data
      for ( int y=0; y < dims[1]; ++y )
	result->Set( pixelBufferTo[y], ThisConst->m_DataGrid->GetOffsetFromIndex( x, y, z ) );
      }
    }
}

void
DataGridFilter
::GetFilteredDataThreadZ( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{ 
  FilterThreadParameters* params = static_cast<FilterThreadParameters*>( args );
  const Self* ThisConst = params->thisObject;
  
  const DataGrid::IndexType& dims = ThisConst->m_DataGrid->m_Dims;
  unsigned int maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const size_t filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray::SmartPtr& result = params->m_Result;

  for ( int y=taskIdx; y < dims[1]; y+=taskCnt ) 
    {
    for ( int x=0; x < dims[0]; ++x ) 
      {
      for ( int z=0; z < dims[2]; ++z )
	if ( !result->Get( pixelBufferFrom[z], ThisConst->m_DataGrid->GetOffsetFromIndex( x, y, z ) ) ) 
	  pixelBufferFrom[z] = 0;
      
      for ( int z=0; z < dims[2]; ++z ) 
	{
	Types::DataItem correctOverlap = filter[0];
	pixelBufferTo[z] = pixelBufferFrom[z] * filter[0];
	for ( int t=1; t < (int)filterSize; ++t ) 
	  {
	  if ( z >= t )
	    {
	    pixelBufferTo[z] += pixelBufferFrom[z-t] * filter[t];
	    correctOverlap += filter[t];
	    }
	  
	  if ( z+t < dims[2] )
	    {
	    pixelBufferTo[z] += pixelBufferFrom[z+t] * filter[t];
	    correctOverlap += filter[t];
	    }
	  }
	// correct value scaling for all missing (outside) elements
	pixelBufferTo[z] /= correctOverlap;
	}
      
      // write back convolved data
      for ( int z=0; z < dims[2]; ++z )
	result->Set( pixelBufferTo[z], ThisConst->m_DataGrid->GetOffsetFromIndex( x, y, z ) );
      }
    }
}

} // namespace cmtk
