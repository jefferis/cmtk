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

#include <cmtkDataGrid.h>
#include <cmtkMemory.h>
#include <cmtkThreadPool.h>

#include <vector>

namespace
cmtk
{

/** \addtogroup Base */
//@{

TypedArray* 
DataGrid::GetFilteredData
( const std::vector<Types::DataItem>& filterX, 
  const std::vector<Types::DataItem>& filterY,
  const std::vector<Types::DataItem>& filterZ ) const
{
  if ( ! this->m_Data ) return NULL;

  TypedArray *result = this->m_Data->NewTemplateArray( this->m_Data->GetDataSize() );

  const size_t numberOfTasks = 4 * ThreadPool::GlobalThreadPool.GetNumberOfThreads() - 3;
  std::vector<FilterThreadParameters> params( numberOfTasks );
  
  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    params[task].thisObject = this;
    params[task].m_Filter = &filterX;
    params[task].m_Result = result;
    }
  ThreadPool::GlobalThreadPool.Run( GetFilteredDataThreadX, params );

  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    params[task].m_Filter = &filterY;
    }
  ThreadPool::GlobalThreadPool.Run( GetFilteredDataThreadY, params );

  for ( size_t task = 0; task < numberOfTasks; ++task )
    {
    params[task].m_Filter = &filterZ;
    }
  ThreadPool::GlobalThreadPool.Run( GetFilteredDataThreadZ, params );
  
  return result;
}

void
DataGrid
::GetFilteredDataThreadX( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  FilterThreadParameters* params = static_cast<FilterThreadParameters*>( args );
  const Self* ThisConst = params->thisObject;
  
  const int* dims = ThisConst->m_Dims;
  unsigned int maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const size_t filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray* result = params->m_Result;

  for ( int z=taskIdx; z < dims[2]; z += taskCnt ) 
    {
    for ( int y=0; y < dims[1]; ++y ) 
      {
      // copy row data to buffer
      size_t ofs = ThisConst->GetOffsetFromIndex( 0, y, z );
      for ( int x=0; x < dims[0]; ++x, ++ofs )
	if ( !ThisConst->m_Data->Get( pixelBufferFrom[x], ofs ) )
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
      
      ofs = ThisConst->GetOffsetFromIndex( 0, y, z );
      for ( int x=0; x < dims[0]; ++x, ++ofs )
	result->Set( pixelBufferTo[x], ofs );
      }
    }
}

void
DataGrid
::GetFilteredDataThreadY( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  FilterThreadParameters* params = static_cast<FilterThreadParameters*>( args );
  const Self* ThisConst = params->thisObject;
  
  const int* dims = ThisConst->m_Dims;
  unsigned int maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const size_t filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray* result = params->m_Result;

  for ( int z=taskIdx; z < dims[2]; z+=taskCnt ) 
    {
    for ( int x=0; x < dims[0]; ++x ) 
      {
      size_t ofs = ThisConst->GetOffsetFromIndex( x, 0, z );
      for ( int y=0; y < dims[1]; ++y, ofs+=ThisConst->nextJ )
	if ( !result->Get( pixelBufferFrom[y], ofs ) ) 
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
      ofs = ThisConst->GetOffsetFromIndex( x, 0, z );
      for ( int y=0; y < dims[1]; ++y, ofs += ThisConst->nextJ )
	result->Set( pixelBufferTo[y], ofs );
      }
    }
}

void
DataGrid
::GetFilteredDataThreadZ( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{ 
  FilterThreadParameters* params = static_cast<FilterThreadParameters*>( args );
  const Self* ThisConst = params->thisObject;
  
  const int* dims = ThisConst->m_Dims;
  unsigned int maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const size_t filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray* result = params->m_Result;

  for ( int y=taskIdx; y < dims[1]; y+=taskCnt ) 
    {
    for ( int x=0; x < dims[0]; ++x ) 
      {
      size_t ofs = ThisConst->GetOffsetFromIndex( x, y, 0 );
      for ( int z=0; z < dims[2]; ++z, ofs+=ThisConst->nextK )
	if ( !result->Get( pixelBufferFrom[z], ofs ) ) 
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
      ofs = ThisConst->GetOffsetFromIndex( x, y, 0 );
      for ( int z=0; z < dims[2]; ++z, ofs += ThisConst->nextK )
	result->Set( pixelBufferTo[z], ofs );
      }
    }
}

} // namespace cmtk
