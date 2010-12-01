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

#include <Base/cmtkTransformedVolumeAxes.h>

namespace
cmtk
{

template<class VM>
typename ImagePairAffineRegistrationFunctionalTemplate<VM>::ReturnType 
ImagePairAffineRegistrationFunctionalTemplate<VM>
::Evaluate() 
{
  const TransformedVolumeAxes axesHash( *this->m_ReferenceGrid, this->m_AffineXform, this->m_FloatingGrid->Deltas().begin(), this->m_FloatingGrid->m_Offset.begin() );
  const Vector3D *axesHashX = axesHash[0], *axesHashY = axesHash[1], *axesHashZ = axesHash[2];
  
  this->m_Metric->Reset();
  
  const DataGrid::IndexType& Dims = this->m_ReferenceGrid->GetDims();
  const int DimsX = Dims[0], DimsY = Dims[1], DimsZ = Dims[2];
  
  this->Clipper.SetDeltaX( axesHashX[DimsX-1] - axesHashX[0] );
  this->Clipper.SetDeltaY( axesHashY[DimsY-1] - axesHashY[0] );
  this->Clipper.SetDeltaZ( axesHashZ[DimsZ-1] - axesHashZ[0] );
  this->Clipper.SetClippingBoundaries( this->m_FloatingCropRegionFractIndex );
  
  DataGrid::IndexType::ValueType startZ, endZ;
  if ( this->ClipZ( this->Clipper, axesHashZ[0], startZ, endZ ) ) 
    {
    startZ = std::max<DataGrid::IndexType::ValueType>( startZ, this->m_ReferenceCropRegion.From()[2] );
    endZ = std::min<DataGrid::IndexType::ValueType>( endZ, this->m_ReferenceCropRegion.To()[2] + 1 );
    
    const int numberOfTasks = std::min<size_t>( 4 * this->m_NumberOfThreads - 3, endZ - startZ + 1 );
    this->m_EvaluateTaskInfo.resize( numberOfTasks );
    
    for ( int taskIdx = 0; taskIdx < numberOfTasks; ++taskIdx ) 
      {
      this->m_EvaluateTaskInfo[taskIdx].thisObject = this;
      this->m_EvaluateTaskInfo[taskIdx].AxesHash = &axesHash;
      this->m_EvaluateTaskInfo[taskIdx].StartZ = startZ;
      this->m_EvaluateTaskInfo[taskIdx].EndZ = endZ;
      }
    
    ThreadPool::GetGlobalThreadPool().Run( EvaluateThread, this->m_EvaluateTaskInfo );
    }
  return this->m_Metric->Get();
}

template<class VM>
void 
ImagePairAffineRegistrationFunctionalTemplate<VM>
::EvaluateThread( void *const args, const size_t taskIdx, const size_t taskCnt, const size_t threadIdx, const size_t ) 
{
  typename Self::EvaluateTaskInfo *info = static_cast<typename Self::EvaluateTaskInfo*>( args );
  
  Self *me = info->thisObject;
  VM& metric = dynamic_cast<VM&>( *(me->m_Metric) );
  
  VM& threadMetric = me->m_ThreadMetric[threadIdx];
  threadMetric.Reset();
  
  const Vector3D *hashX = (*info->AxesHash)[0], *hashY = (*info->AxesHash)[1], *hashZ = (*info->AxesHash)[2];
  Vector3D pFloating;
  
  const DataGrid::IndexType& Dims = me->m_ReferenceGrid->GetDims();
  const int DimsX = Dims[0], DimsY = Dims[1];
  
  int fltIdx[3];
  Types::Coordinate fltFrac[3];
  
  Vector3D rowStart;
  Vector3D planeStart;
  
  DataGrid::IndexType::ValueType pX, pY, pZ;
  // Loop over all remaining planes
  for ( pZ = info->StartZ + taskIdx; pZ < info->EndZ; pZ += taskCnt ) 
    {
    // Offset of current reference voxel
    int r = pZ * DimsX * DimsY;
    
    planeStart = hashZ[pZ];
    
    DataGrid::IndexType::ValueType startY, endY;
    if ( me->ClipY( me->Clipper, planeStart, startY, endY ) ) 
      {	
      startY = std::max<DataGrid::IndexType::ValueType>( startY, me->m_ReferenceCropRegion.From()[1] );
      endY = std::min<DataGrid::IndexType::ValueType>( endY, me->m_ReferenceCropRegion.To()[1] + 1 );
      r += startY * DimsX;
      
      // Loop over all remaining rows
      for ( pY = startY; pY<endY; ++pY ) 
	{
	(rowStart = planeStart) += hashY[pY];
	
	DataGrid::IndexType::ValueType startX, endX;
	if ( me->ClipX( me->Clipper, rowStart, startX, endX ) ) 
	  {	    
	  startX = std::max<DataGrid::IndexType::ValueType>( startX, me->m_ReferenceCropRegion.From()[0] );
	  endX = std::min<DataGrid::IndexType::ValueType>( endX, me->m_ReferenceCropRegion.To()[0] + 1 );
	  
	  r += startX;
	  // Loop over all remaining voxels in current row
	  for ( pX = startX; pX<endX; ++pX, ++r ) 
	    {
	    // Continue metric computation.
	    Types::DataItem sampleX;
	    if ( metric.GetSampleX( sampleX, r ) )
	      {
	      (pFloating = rowStart) += hashX[pX];
	      
	      // probe volume and get the respective voxel
	      if ( me->m_FloatingGrid->FindVoxelByIndex( pFloating, fltIdx, fltFrac ) )
		{
		threadMetric.Increment( sampleX, metric.GetSampleY( fltIdx, fltFrac ) );
		}
	      else
		{ 
		if ( me->m_ForceOutsideFlag ) 
		  {
		  threadMetric.Increment( sampleX, me->m_ForceOutsideValueRescaled );		
		  }
		}
	      }
	    }
	  r += (DimsX-endX);
	  } 
	else
	  {
	  r += DimsX;
	  }
	}
      
      r += (DimsY-endY) * DimsX;
      } 
    else
      {
      r += DimsY * DimsX;
      }
    }
  
  me->m_MetricMutex.Lock();
  metric.Add( threadMetric );
  me->m_MetricMutex.Unlock();
}

} // namnespace cmtk
