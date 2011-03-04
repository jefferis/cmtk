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

//#define DEBUG_COMM

#include <Registration/cmtkGroupwiseRegistrationFunctionalXformTemplate.h>

#include <Base/cmtkMathUtil.h>
#include <IO/cmtkVolumeIO.h>

#ifdef CMTK_BUILD_MPI
#  include <mpi.h>
#  include <IO/cmtkMPI.h>
#endif

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TXform>
GroupwiseRegistrationFunctionalXformTemplate<TXform>::GroupwiseRegistrationFunctionalXformTemplate()
{
  // allocate storage for four times me parallel tasks than threads in the global thread pool.
  this->m_InterpolateTaskInfo.resize( 4 * ThreadPool::GetGlobalThreadPool().GetNumberOfThreads() );
}

template<class TXform>
void
GroupwiseRegistrationFunctionalXformTemplate<TXform>
::InterpolateImage
( const size_t idx, byte* const destination )
{ 
  for ( size_t task = 0; task < this->m_InterpolateTaskInfo.size(); ++task )
    {
    this->m_InterpolateTaskInfo[task].thisObject = this;
    this->m_InterpolateTaskInfo[task].m_Idx = idx;    
    this->m_InterpolateTaskInfo[task].m_Destination = destination;    
    }

  ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
  if ( this->m_ProbabilisticSamples.size() )
    threadPool.Run( InterpolateImageProbabilisticThread, this->m_InterpolateTaskInfo );
  else
    threadPool.Run( InterpolateImageThread, this->m_InterpolateTaskInfo );

  // Sum number of pixels outside FOV from all tasks.
  size_t numberOfOutsidePixels = 0;
  for ( size_t task = 0; task < this->m_InterpolateTaskInfo.size(); ++task )
    {
    numberOfOutsidePixels += this->m_InterpolateTaskInfo[task].m_NumberOfOutsidePixels;
    }

  // Check whether more than defined proportion threshold was outside
  if ( numberOfOutsidePixels > static_cast<size_t>( this->m_MaxRelativeNumberOutsidePixels * this->m_TemplateNumberOfSamples ) )
    {
    throw typename Self::BadXform();
    }
}

template<class TXform>
void
GroupwiseRegistrationFunctionalXformTemplate<TXform>::InterpolateImageThread
( void *const args, const size_t taskIdx, const size_t taskCnt )
{
  InterpolateImageThreadParameters* threadParameters = static_cast<InterpolateImageThreadParameters*>( args );
  
  const Self* This = threadParameters->thisObject;
  const size_t idx = threadParameters->m_Idx;
  byte* destination = threadParameters->m_Destination;

  const TXform* xform = This->GetXformByIndex(idx);
  const UniformVolume* target = This->m_ImageVector[idx];

  const byte paddingValue = This->m_PaddingValue;
  const byte backgroundValue = This->m_UserBackgroundFlag ? This->m_PrivateUserBackgroundValue : paddingValue;
  threadParameters->m_NumberOfOutsidePixels = 0;

  Vector3D v;
  byte value;
  const byte* dataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );

  const int dimsX = This->m_TemplateGrid->GetDims()[AXIS_X];
  const int dimsY = This->m_TemplateGrid->GetDims()[AXIS_Y];
  const int dimsZ = This->m_TemplateGrid->GetDims()[AXIS_Z];

  for ( int z = taskIdx; (z < dimsZ); z += taskCnt ) 
    {
    byte *wptr = destination + z * dimsX * dimsY;
    for ( int y = 0; (y < dimsY); ++y )
      {
      for ( int x = 0; x < dimsX; ++x )
	{
	This->m_TemplateGrid->GetGridLocation( v, x, y, z );
	xform->ApplyInPlace( v );
	
	if ( target->ProbeData( value, dataPtr, v ) )
	  {
	  *wptr = value;
	  }
	else
	  {
	  *wptr = backgroundValue;
	  ++threadParameters->m_NumberOfOutsidePixels;
	  }
	
	++wptr;
	}
      }
    }
}

template<class TXform>
void
GroupwiseRegistrationFunctionalXformTemplate<TXform>::InterpolateImageProbabilisticThread
( void *const args, const size_t taskIdx, const size_t taskCnt )
{
  InterpolateImageThreadParameters* threadParameters = static_cast<InterpolateImageThreadParameters*>( args );
  
  const Self* This = threadParameters->thisObject;
  const size_t idx = threadParameters->m_Idx;
  byte* destination = threadParameters->m_Destination;
  threadParameters->m_NumberOfOutsidePixels = 0;

  const TXform* xform = This->GetXformByIndex(idx);
  const UniformVolume* target = This->m_ImageVector[idx];

  const byte paddingValue = This->m_PaddingValue;
  const byte backgroundValue = This->m_UserBackgroundFlag ? This->m_PrivateUserBackgroundValue : paddingValue;

  Vector3D v;
  byte value;
  const byte* dataPtr = static_cast<const byte*>( target->GetData()->GetDataPtr() );

  const size_t startIdx = taskIdx * (This->m_ProbabilisticSamples.size() / taskCnt);
  const size_t endIdx = ( taskIdx == taskCnt ) ? This->m_ProbabilisticSamples.size() : (taskIdx+1) * (This->m_ProbabilisticSamples.size() / taskCnt);

  byte *wptr = destination + startIdx;
  for ( size_t i = startIdx; i < endIdx; ++i, ++wptr )
    {
    const size_t offset = This->m_ProbabilisticSamples[i];
    This->m_TemplateGrid->GetGridLocation( v, offset );
    xform->ApplyInPlace( v );
    
    if ( target->ProbeData( value, dataPtr, v ) )
      {
      *wptr = value;
      }
    else
      {
      *wptr = backgroundValue;
      ++threadParameters->m_NumberOfOutsidePixels;
      }
    }
}

//@}

} // namespace cmtk
