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

#include <cmtkUniformVolume.h>

#include <cmtkVolumeGridToGridLookup.h>

namespace
cmtk
{

/** \addtogroup Base */
//@{

TypedArray*
UniformVolume::Resample( const UniformVolume& other ) const 
{
  const TypedArray* fromData = other.GetData();

  const VolumeGridToGridLookup gridLookup( other, *this );
  const size_t numberOfThreads = std::min<int>( Threads::GetNumberOfThreads(), this->m_Dims[2] );
   
  // Info blocks for parallel threads that do the resampling.
  UniformVolume::ResampleThreadInfo *ThreadInfo = Memory::AllocateArray<UniformVolume::ResampleThreadInfo>( numberOfThreads );
   
  Types::DataItem *resampledData = NULL;
  try
    {
    resampledData = Memory::AllocateArray<Types::DataItem>(  this->GetNumberOfPixels()  );
     
    for ( size_t threadIdx = 0; threadIdx < numberOfThreads; ++threadIdx ) 
      {
          std::cerr << threadIdx << "\t" << &ThreadInfo[threadIdx] << std::endl;
          std::cerr.flush();
      ThreadInfo[threadIdx].thisObject = this;
      ThreadInfo[threadIdx].ThisThreadIndex = threadIdx;
      ThreadInfo[threadIdx].NumberOfThreads = numberOfThreads;
      ThreadInfo[threadIdx].GridLookup = &gridLookup;
      ThreadInfo[threadIdx].OtherVolume = &other;
      ThreadInfo[threadIdx].FromData = fromData;
      ThreadInfo[threadIdx].ResampledData = resampledData;
      }
     
    switch ( fromData->GetDataClass() ) 
      {
      case DATACLASS_BINARY:
	fputs( "Data class 'binary' not supported in UniformVolume::Resample.\n"
	       "Falling back to 'greyscale'. "
	       "You may not like the result though...\n", stderr );
      case DATACLASS_GREY:
      default:
	Threads::RunThreads( UniformVolume::ResampleThreadExecuteGrey, numberOfThreads, ThreadInfo );
	break;
      case DATACLASS_LABEL:
	Threads::RunThreads( UniformVolume::ResampleThreadExecuteLabels, numberOfThreads, ThreadInfo );
	break;
      }
    }
   
  catch ( std::bad_alloc ) 
    {
    resampledData = NULL;
    }
   
  delete[] ThreadInfo;
   
  if ( resampledData == NULL ) 
    {
    return NULL;
    }
   
  TypedArray *Result = fromData->NewTemplateArray( resampledData, this->GetNumberOfPixels() );
  delete[] resampledData;
   
  Result->SetDataClass( fromData->GetDataClass() );
   
  return Result;
}

CMTK_THREAD_RETURN_TYPE
UniformVolume::ResampleThreadExecuteLabels( void *arg )
{
  UniformVolume::ResampleThreadInfo *info = static_cast<UniformVolume::ResampleThreadInfo*>( arg );

  const UniformVolume *me = info->thisObject;
  const UniformVolume *other = info->OtherVolume;
  Types::DataItem *dest = info->ResampledData;
  const VolumeGridToGridLookup *gridLookup = info->GridLookup;
  
  Types::Coordinate labelWeights[256];

  int x, y;
  int pX, pY, pZ;
  Types::DataItem value;
  
  for ( int z = info->ThisThreadIndex; z < me->m_Dims[2]; z += info->NumberOfThreads ) 
    {
    int offset = z * me->m_Dims[0] * me->m_Dims[1];
    
    for ( y = 0; y < me->m_Dims[1]; ++y ) 
      {
      for ( x = 0; x < me->m_Dims[0]; ++x, ++offset ) 
	{
	memset( labelWeights, 0, sizeof( labelWeights ) );
	
	for ( pZ=0; pZ<gridLookup->GetSourceCount(2,z); ++pZ ) 
	  {
		  const Types::Coordinate weightZ=gridLookup->GetWeight(2,z,pZ);
	  
	  for ( pY=0; pY<gridLookup->GetSourceCount(1,y); ++pY ) 
	    {
	    const Types::Coordinate weightYZ=weightZ*gridLookup->GetWeight(1,y,pY);
	    
	    for ( pX=0; pX<gridLookup->GetSourceCount(0,x); ++pX ) 
	      {
	      const Types::Coordinate weight=weightYZ*gridLookup->GetWeight(0,x,pX);
	      
	      if ( other->GetDataAt( value, pX + gridLookup->GetFromIndex(0,x), pY + gridLookup->GetFromIndex(1,y), pZ + gridLookup->GetFromIndex(2,z)) )
		{
		labelWeights[static_cast<byte>( value )] += weight;
		} 
	      }
	    }
	  }
	
	Types::Coordinate maxLabelWeight = 0;
	byte maxLabelIndex = 0;
	for ( int l=0; l<256; ++l ) 
	  {
	  if ( labelWeights[l] > maxLabelWeight ) 
	    {
	    maxLabelWeight = labelWeights[l];
	    maxLabelIndex = l;
	    }
	  }
	
	if ( maxLabelWeight > 0 )
	  dest[offset] = maxLabelIndex;
	else 
	  dest[offset] = sqrt( static_cast<Types::DataItem>( -1 ) );
	}
      }
    }
  return CMTK_THREAD_RETURN_VALUE;
}

CMTK_THREAD_RETURN_TYPE
UniformVolume::ResampleThreadExecuteGrey( void *arg )
{
  UniformVolume::ResampleThreadInfo *info = static_cast<UniformVolume::ResampleThreadInfo*>( arg );

  const UniformVolume *me = info->thisObject;
  Types::DataItem *dest = info->ResampledData;
  const UniformVolume *other = info->OtherVolume;
  const VolumeGridToGridLookup *gridLookup = info->GridLookup;
  
  Types::DataItem tempValue, value;
  
  int x, y;
  int pX, pY, pZ;
  bool FoundNullData;
  
  for ( int z = info->ThisThreadIndex; z < me->m_Dims[2]; z += info->NumberOfThreads ) 
    {
    int offset = z * me->m_Dims[0] * me->m_Dims[1];

    const Types::Coordinate volumeZ = gridLookup->GetLength(2,z);
	
    for ( y=0; y < me->m_Dims[1]; ++y ) 
      {
      const Types::Coordinate volumeYZ = volumeZ * gridLookup->GetLength(1,y);
      
      for ( x=0; x < me->m_Dims[0]; ++x, ++offset ) 
	{
	tempValue = 0;
	FoundNullData = false;
	
	for ( pZ=0; pZ<gridLookup->GetSourceCount(2,z); ++pZ ) 
	  {
	  const Types::Coordinate weightZ=gridLookup->GetWeight(2,z,pZ);
	  
	  for ( pY=0; pY<gridLookup->GetSourceCount(1,y); ++pY ) 
	    {
	    const Types::Coordinate weightYZ=weightZ*gridLookup->GetWeight(1,y,pY);
	    
	    for ( pX=0; pX<gridLookup->GetSourceCount(0,x); ++pX ) 
	      {
	      const Types::Coordinate weight=weightYZ*gridLookup->GetWeight(0,x,pX);
	      
	      if ( other->GetDataAt(value,pX+gridLookup->GetFromIndex(0,x), pY+gridLookup->GetFromIndex(1,y), pZ+gridLookup->GetFromIndex(2,z) ) )
		{
			tempValue+=static_cast<Types::DataItem>( weight*value );
		} 
	      else
		{
		FoundNullData = true;
		}
	      }
	    }
	  }
	
	if ( ! FoundNullData ) 
	  {
	  const Types::Coordinate volume = volumeYZ*gridLookup->GetLength(0,x);
	  dest[offset] = static_cast<Types::DataItem>( tempValue / volume );
	  } 
	else
	  {
	  dest[offset] = sqrt( static_cast<Types::DataItem>( -1 ) );
	  }
	}
      }
    }
  return CMTK_THREAD_RETURN_VALUE;
}

} // namespace cmtk
