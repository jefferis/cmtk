/*
//
//  Copyright 2016 Google, Inc.
//
//  Copyright 1997-2009 Torsten Rohlfing
//
//  Copyright 2004-2013 SRI International
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

#include "cmtkDataGridFilter.h"

#include <System/cmtkException.h>
#include <System/cmtkThreadPool.h>
#include <System/cmtkProgress.h>

#include <Base/cmtkRegionIndexIterator.h>
#include <Base/cmtkMathFunctionWrappers.h>

namespace
cmtk
{

DataGridFilter
::DataGridFilter( DataGrid::SmartConstPtr dataGrid ) : m_DataGrid( dataGrid ) {}

cmtk::Types::DataItem
cmtk::DataGridFilter::MedianOperator::Reduce( std::vector<Types::DataItem>& regionData )
{
  std::sort( regionData.begin(), regionData.end() );
  
  if ( regionData.size() % 2 )
    return regionData[regionData.size()/2];
  else
    return (Types::DataItem) ( 0.5 * (regionData[regionData.size()/2] + regionData[regionData.size()/2-1]) );
}

TypedArray::SmartPtr
DataGridFilter::GetDataMedianFiltered( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  return this->ApplyRegionFilter<Self::MedianOperator>( radiusX, radiusY, radiusZ );
}

cmtk::Types::DataItem
cmtk::DataGridFilter::MeanOperator::Reduce( std::vector<Types::DataItem>& regionData )
{
  Types::DataItem sum = 0;
  for ( Types::GridIndexType i = 0; i < regionData.size(); ++i )
    sum += regionData[i];
  
  return sum / regionData.size();
}

TypedArray::SmartPtr
DataGridFilter::RegionMeanFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  return this->ApplyRegionFilter<Self::MeanOperator>( radiusX, radiusY, radiusZ );
}

TypedArray::SmartPtr
DataGridFilter::FastRegionMeanFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  const DataGrid& dataGrid = *(this->m_DataGrid);
  if ( !dataGrid.GetData() )
    return TypedArray::SmartPtr( NULL );

  const TypedArray& dataArray = *(dataGrid.GetData());  

  DataGrid::IndexType radius;
  radius[0] = radiusX;
  radius[1] = radiusY;
  radius[2] = radiusZ;

  const Types::GridIndexType nPixels = dataGrid.GetNumberOfPixels();
  const DataGrid::RegionType wholeImageRegion = dataGrid.GetWholeImageRegion();
  
  std::vector<double> sums( nPixels );
  std::fill( sums.begin(), sums.end(), 0 );

  std::vector<unsigned short> cnts( nPixels );
  std::fill( cnts.begin(), cnts.end(), 0 );

  //
  // Mean filter is separable - process x,y,z separately
  //
  for ( int dim = 0; dim < 3; ++dim )
    {
    const DataGrid::RegionType face = wholeImageRegion.GetFaceRegion( dim );

    const Types::GridIndexType columnFrom = wholeImageRegion.From()[dim];
    const Types::GridIndexType columnTo = wholeImageRegion.To()[dim];
    const Types::GridIndexType nPixelsColumn = columnTo - columnFrom;

    std::vector<double> sumsColumn( nPixelsColumn );
    std::vector<unsigned short> cntsColumn( nPixelsColumn );
    
    for ( RegionIndexIterator<DataGrid::RegionType> fIt = RegionIndexIterator<DataGrid::RegionType>( face ); fIt != fIt.end(); ++fIt )
      {
      double sum = 0;
      unsigned short count = 0;

      //
      // PASS 1 - accumulate all values and counts upwards
      //
      Types::GridIndexType idx0 = 0;
      DataGrid::IndexType idx = fIt.Index();
      
      for ( idx[dim] = columnFrom; idx[dim] < columnTo; ++idx[dim], ++idx0 )
	{
	const Types::GridIndexType offset = dataGrid.GetOffsetFromIndex( idx );

	if ( ! dim ) // get first-pass values from data array
	  {
	  Types::DataItem value;
	  if ( dataArray.Get( value, offset ) )
	    {
	    ++count;
	    }
	  else
	    {
	    value = 0;
	    }
	  cntsColumn[idx0] = count;
	  sumsColumn[idx0] = (sum+=value);
	  }
	else
	  {
	  cntsColumn[idx0] = (count+=cnts[offset]);
	  sumsColumn[idx0] = (sum+=sums[offset]);
	  }
	}
      
      //
      // PASS 2 - compute differences between upper and lower end of kernel window
      //
      idx0 = 0;
      for ( idx[dim] = columnFrom; idx[dim] < columnTo; ++idx[dim], ++idx0 )
	{
	const Types::GridIndexType offset = dataGrid.GetOffsetFromIndex( idx );
	
	const Types::GridIndexType upper = static_cast<Types::GridIndexType>( std::min<Types::GridIndexType>( idx0+radius[dim], nPixelsColumn-1 ) );	
	sums[offset] = sumsColumn[upper];
	cnts[offset] = cntsColumn[upper];

	// if lower end of window is outside range, implicitly subtract zero from sums and counts
	const Types::GridIndexType lower = idx0-radius[dim]-1;
	if ( lower >= 0 )
	  {
	  sums[offset] -= sumsColumn[lower];
	  cnts[offset] -= cntsColumn[lower];
	  } 
	}
      }
    }

  TypedArray::SmartPtr result( TypedArray::Create( dataArray.GetType(), nPixels ) );
  for ( Types::GridIndexType ofs = 0; ofs < nPixels; ++ofs )
    {
    if ( cnts[ofs] )
      {
      result->Set( sums[ofs] / cnts[ofs], ofs );
      }
    else
      {
      result->SetPaddingAt( ofs );
      }
    }
  return result;
}

TypedArray::SmartPtr
DataGridFilter::FastRegionVarianceFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  if ( !this->m_DataGrid->GetData() )
    return TypedArray::SmartPtr( NULL );

  //
  // see http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  //
  TypedArray::SmartPtr mean = this->FastRegionMeanFilter( radiusX, radiusY, radiusZ );

  DataGrid::SmartPtr square( this->m_DataGrid->Clone() );
  square->GetData()->ApplyFunctionDouble( Wrappers::Square );
  square->SetData( Self( square ).FastRegionMeanFilter( radiusX, radiusY, radiusZ ) );

  TypedArray& squareData = *(square->GetData());  

  const Types::GridIndexType nPixels = square->GetNumberOfPixels();  
  for ( Types::GridIndexType i = 0; i < nPixels; ++i )
    {
    Types::DataItem vMean, vSquare;
    if ( mean->Get( vMean, i ) && squareData.Get( vSquare, i ) )
      {
      squareData.Set( vSquare - vMean * vMean, i );
      }
    else
      {
      squareData.SetPaddingAt( i );
      }
    }

  return square->GetData();
}


cmtk::Types::DataItem
cmtk::DataGridFilter::VarianceOperator::Reduce( std::vector<Types::DataItem>& regionData )
{
  const Types::DataItem mean = MeanOperator::Reduce( regionData );

  Types::DataItem sum = 0;
  for ( Types::GridIndexType i = 0; i < regionData.size(); ++i )
    {
    const Types::DataItem diff = mean - regionData[i];
    sum += diff*diff;
    }
  
  return sum / regionData.size();
}

TypedArray::SmartPtr
DataGridFilter::RegionVarianceFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  return this->ApplyRegionFilter<Self::VarianceOperator>( radiusX, radiusY, radiusZ );
}

cmtk::Types::DataItem
cmtk::DataGridFilter::StandardDeviationOperator::Reduce( std::vector<Types::DataItem>& regionData )
{
  return sqrt( VarianceOperator::Reduce( regionData ) );
}

TypedArray::SmartPtr
DataGridFilter::RegionStandardDeviationFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  return this->ApplyRegionFilter<Self::StandardDeviationOperator>( radiusX, radiusY, radiusZ );
}

cmtk::Types::DataItem
cmtk::DataGridFilter::SmoothnessOperator::Reduce( std::vector<Types::DataItem>& regionData )
{
  return (1.0 - 1.0/(1.0 + VarianceOperator::Reduce( regionData ) ));
}

TypedArray::SmartPtr
DataGridFilter::RegionSmoothnessFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  return this->ApplyRegionFilter<Self::SmoothnessOperator>( radiusX, radiusY, radiusZ );
}

cmtk::Types::DataItem
cmtk::DataGridFilter::ThirdMomentOperator::Reduce( std::vector<Types::DataItem>& regionData )
{
  const Types::DataItem mean = MeanOperator::Reduce( regionData );

  Types::DataItem sum = 0.0;
  for ( Types::GridIndexType i = 0; i < regionData.size(); ++i ) 
    {
    const Types::DataItem diff = mean - regionData[i];
    sum += diff*diff*diff;
    }

  return sum / MathUtil::Square( regionData.size() );
}

TypedArray::SmartPtr
DataGridFilter::RegionThirdMomentFilter( const Types::GridIndexType radiusX, const Types::GridIndexType radiusY, const Types::GridIndexType radiusZ ) const
{
  return this->ApplyRegionFilter<Self::ThirdMomentOperator>( radiusX, radiusY, radiusZ );
}

TypedArray::SmartPtr
DataGridFilter::GetDataSobelFiltered() const
{
  const TypedArray* data = this->m_DataGrid->GetData();
  if ( !data )
    return TypedArray::SmartPtr( NULL );

  TypedArray::SmartPtr result = TypedArray::Create( data->GetType(), data->GetDataSize() );

  Types::DataItem value = 0;
  Types::DataItem fov[3][3][3];

  Progress::Begin( 0, this->m_DataGrid->m_Dims[2], 1 );

  Types::GridIndexType offset = 0;
  for ( Types::GridIndexType z = 0; z < this->m_DataGrid->m_Dims[2]; ++z ) 
    {
    Progress::SetProgress( z );
    for ( Types::GridIndexType y = 0; y < this->m_DataGrid->m_Dims[1]; ++y )
      for ( Types::GridIndexType x = 0; x < this->m_DataGrid->m_Dims[0]; ++x, ++offset ) 
	{
	if ( x && y && z && ( x<this->m_DataGrid->m_Dims[0]-1 ) && ( y<this->m_DataGrid->m_Dims[1]-1 ) && ( z<this->m_DataGrid->m_Dims[2]-1 ) ) 
	  {
	  for ( Types::GridIndexType dz=-1; dz<2; ++dz )
	    for ( Types::GridIndexType dy=-1; dy<2; ++dy )
	      for ( Types::GridIndexType dx=-1; dx<2; ++dx )
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
  const std::vector<Types::DataItem>& filterZ,
  const bool normalize ) const
{
  if ( !this->m_DataGrid->GetData() )
    return TypedArray::SmartPtr( NULL );

  // create and initialize result (data are copied at this point from the source array)
  TypedArray::SmartPtr result( this->m_DataGrid->GetData()->Clone() );
  
  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  const Types::GridIndexType numberOfTasks = 4 * threadPool.GetNumberOfThreads() - 3;

  // initialize common thread data
  std::vector<FilterThreadParameters> params( numberOfTasks );
  for ( Types::GridIndexType task = 0; task < numberOfTasks; ++task )
    {
    params[task].thisObject = this;
    params[task].m_Result = result;
    params[task].m_Normalize = normalize;
    }

  // run x filter
  if ( filterX.size() > 1 )
    {
    for ( Types::GridIndexType task = 0; task < numberOfTasks; ++task )
      {
      params[task].m_Filter = &filterX;
      }
    threadPool.Run( GetFilteredDataThreadX, params );
    }

  // run y filter
  if ( filterY.size() > 1 )
    {
    for ( Types::GridIndexType task = 0; task < numberOfTasks; ++task )
      {
      params[task].m_Filter = &filterY;
      }
    threadPool.Run( GetFilteredDataThreadY, params );
    }
  
  // run z filter
  if ( filterZ.size() > 1 )
    {
    for ( Types::GridIndexType task = 0; task < numberOfTasks; ++task )
      {
      params[task].m_Filter = &filterZ;
      }
    threadPool.Run( GetFilteredDataThreadZ, params );
    }
  
  return result;
}

void
DataGridFilter
::GetFilteredDataThreadX( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t )
{
  FilterThreadParameters* params = static_cast<FilterThreadParameters*>( args );
  const Self* ThisConst = params->thisObject;
  
  const DataGrid::IndexType& dims = ThisConst->m_DataGrid->m_Dims;
  Types::GridIndexType maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const bool normalize = params->m_Normalize;
  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const Types::GridIndexType filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray::SmartPtr& result = params->m_Result;

  for ( Types::GridIndexType z=taskIdx; z < dims[2]; z += taskCnt ) 
    {
    for ( Types::GridIndexType y=0; y < dims[1]; ++y ) 
      {
      // copy row data to buffer
      Types::GridIndexType ofs = ThisConst->m_DataGrid->GetOffsetFromIndex( 0, y, z );
      for ( Types::GridIndexType x=0; x < dims[0]; ++x, ++ofs )
	if ( !result->Get( pixelBufferFrom[x], ofs ) )
	  pixelBufferFrom[x] = 0;
      
      // convolve row with filter
      for ( Types::GridIndexType x=0; x < dims[0]; ++x ) 
	{
	// this keeps track of outside FOV data
	Types::DataItem correctOverlap = filter[0];
	// central element first to initialize target value
	pixelBufferTo[x] = pixelBufferFrom[x] * filter[0];
	// now convolve side elements
	for ( Types::GridIndexType t=1; t < (Types::GridIndexType)filterSize; ++t ) 
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
	if ( normalize && (correctOverlap != 0) )
	  pixelBufferTo[x] /= correctOverlap;
	}
      
      ofs = ThisConst->m_DataGrid->GetOffsetFromIndex( 0, y, z );
      for ( Types::GridIndexType x=0; x < dims[0]; ++x, ++ofs )
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
  Types::GridIndexType maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const bool normalize = params->m_Normalize;
  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const Types::GridIndexType filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray::SmartPtr& result = params->m_Result;

  for ( Types::GridIndexType z=taskIdx; z < dims[2]; z+=taskCnt ) 
    {
    for ( Types::GridIndexType x=0; x < dims[0]; ++x ) 
      {
      for ( Types::GridIndexType y=0; y < dims[1]; ++y )
	if ( !result->Get( pixelBufferFrom[y], ThisConst->m_DataGrid->GetOffsetFromIndex( x, y, z ) ) ) 
	  pixelBufferFrom[y] = 0;
      
      for ( Types::GridIndexType y=0; y < dims[1]; ++y ) 
	{
	Types::DataItem correctOverlap = filter[0];
	pixelBufferTo[y] = pixelBufferFrom[y] * filter[0];
	for ( Types::GridIndexType t=1; t < (Types::GridIndexType)filterSize; ++t ) 
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
	if ( normalize && (correctOverlap != 0) )
	  pixelBufferTo[y] /= correctOverlap;
	}
      
      // write back convolved data
      for ( Types::GridIndexType y=0; y < dims[1]; ++y )
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
  Types::GridIndexType maxDim = std::max( dims[0], std::max( dims[1], dims[2] ) );

  const bool normalize = params->m_Normalize;
  const std::vector<Types::DataItem>& filter = *(params->m_Filter);
  const Types::GridIndexType filterSize = filter.size();

  std::vector<Types::DataItem> pixelBufferFrom( maxDim );
  std::vector<Types::DataItem> pixelBufferTo( maxDim );
  TypedArray::SmartPtr& result = params->m_Result;

  for ( Types::GridIndexType y=taskIdx; y < dims[1]; y+=taskCnt ) 
    {
    for ( Types::GridIndexType x=0; x < dims[0]; ++x ) 
      {
      for ( Types::GridIndexType z=0; z < dims[2]; ++z )
	if ( !result->Get( pixelBufferFrom[z], ThisConst->m_DataGrid->GetOffsetFromIndex( x, y, z ) ) ) 
	  pixelBufferFrom[z] = 0;
      
      for ( Types::GridIndexType z=0; z < dims[2]; ++z ) 
	{
	Types::DataItem correctOverlap = filter[0];
	pixelBufferTo[z] = pixelBufferFrom[z] * filter[0];
	for ( Types::GridIndexType t=1; t < (Types::GridIndexType)filterSize; ++t ) 
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
	if ( normalize && (correctOverlap != 0) )
	  pixelBufferTo[z] /= correctOverlap;
	}
      
      // write back convolved data
      for ( Types::GridIndexType z=0; z < dims[2]; ++z )
	result->Set( pixelBufferTo[z], ThisConst->m_DataGrid->GetOffsetFromIndex( x, y, z ) );
      }
    }
}

} // namespace cmtk
