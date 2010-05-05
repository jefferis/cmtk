/*
//
//  Copyright 1997-2010 Torsten Rohlfing
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

#include <cmtkThreadPool.h>

namespace
cmtk
{

/** \addtogroup Registration */
//@{

template<class TMultiChannelMetricFunctional>
void
AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional>
::InitTransformation( const bool alignCenters )
{
  this->m_Transformation.MakeIdentityXform();
  if ( alignCenters ) 
    {
    if ( this->m_ReferenceChannels.size() == 0 || this->m_FloatingChannels.size() == 0 )
      {
      StdErr << "ERROR: must set at least one reference and one floating channel image before calling\n"
		<< "       AffineMultiChannelRegistrationFunctional::InitTransformation()\n";
      exit( 1 );
      }

    Vector3D deltaCenter = ( this->m_ReferenceChannels[0]->GetCenterCropRegion() - this->m_FloatingChannels[0]->GetCenterCropRegion() );
    this->m_Transformation.SetXlate( deltaCenter.XYZ );
    }
  
  Vector3D center = this->m_ReferenceChannels[0]->GetCenterCropRegion();
  this->m_Transformation.ChangeCenter( center.XYZ );
}

template<class TMultiChannelMetricFunctional>
typename AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional>::ReturnType
AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional>
::Evaluate() 
{
  this->m_MetricData.Init( this );

  const VolumeAxesHash transformedAxes( *this->m_ReferenceChannels[0], &this->m_Transformation );
  
  const DataGrid::IndexType& dims = this->m_ReferenceDims;
  const int dimsX = dims[0], dimsY = dims[1], dimsZ = dims[2];

  this->m_VolumeClipper.SetDeltaX( transformedAxes[0][dimsX-1] - transformedAxes[0][0] );
  this->m_VolumeClipper.SetDeltaY( transformedAxes[1][dimsY-1] - transformedAxes[1][0] );
  this->m_VolumeClipper.SetDeltaZ( transformedAxes[2][dimsZ-1] - transformedAxes[2][0] );
  this->m_VolumeClipper.SetClippingBoundaries( this->m_FloatingCropRegion );
    
  int startZ, endZ;
  if ( this->ClipZ( this->m_VolumeClipper, transformedAxes[2][0], startZ, endZ ) ) 
    {
    startZ = std::max<int>( startZ, this->m_ReferenceCropRegion.From()[2] );
    endZ = std::min<int>( endZ, this->m_ReferenceCropRegion.To()[2] );

    ThreadPool& threadPool = ThreadPool::GetGlobalThreadPool();
    const size_t numberOfThreads = threadPool.GetNumberOfThreads();
    const size_t numberOfTasks = 4 * numberOfThreads - 3;

    std::vector<typename Self::EvaluateThreadParameters> threadParams( numberOfTasks );
    for ( size_t thread = 0; thread < numberOfTasks; ++thread )
      {      
      threadParams[thread].thisObject = this;
      threadParams[thread].m_TransformedAxes = &transformedAxes;
      threadParams[thread].m_StartZ = startZ;
      threadParams[thread].m_EndZ = endZ;
      }
    threadPool.Run( Self::EvaluateThreadFunction, threadParams );
    }
    
  return this->GetMetric( this->m_MetricData );
}
      
template<class TMultiChannelMetricFunctional>
void
AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional>
::EvaluateThreadFunction( void* args, const size_t taskIdx, const size_t taskCnt, const size_t, const size_t  ) 
{
  typename Self::EvaluateThreadParameters* params = static_cast<typename Self::EvaluateThreadParameters*>( args );
  
  Self* This = params->thisObject;
  const Self* constThis = This;
  
  typename Self::MetricData metricData;
  metricData.Init( This );
  
  const Vector3D *transformedAxesX = (*params->m_TransformedAxes)[0];
  const Vector3D *transformedAxesY = (*params->m_TransformedAxes)[1];
  const Vector3D *transformedAxesZ = (*params->m_TransformedAxes)[2];

  const DataGrid::IndexType& dims = constThis->m_ReferenceDims;
  const int dimsX = dims[0], dimsY = dims[1];

  Vector3D pFloating, rowStart;

  // Loop over all remaining planes
  for ( int pZ = params->m_StartZ + taskIdx; pZ < params->m_EndZ; pZ += taskCnt ) 
    {
    // Offset of current reference voxel
    int r = pZ * dimsX * dimsY;
    Vector3D planeStart = transformedAxesZ[pZ];
    
    int startY, endY;
    if ( constThis->ClipY( constThis->m_VolumeClipper, planeStart, startY, endY ) ) 
      {	  
      startY = std::max<int>( startY, constThis->m_ReferenceCropRegion.From()[1] );
      endY = std::min<int>( endY, constThis->m_ReferenceCropRegion.To()[1] );
      r += startY * dimsX;
      
      // Loop over all remaining rows
      for ( int pY = startY; pY<endY; ++pY ) 
	{
	(rowStart = planeStart) += transformedAxesY[pY];
	
	int startX, endX;
	if ( constThis->ClipX( constThis->m_VolumeClipper, rowStart, startX, endX ) ) 
	  {
	  startX = std::max<int>( startX, constThis->m_ReferenceCropRegion.From()[0] );
	  endX = std::min<int>( endX, constThis->m_ReferenceCropRegion.To()[0] );
	  
	  r += startX;
	  // Loop over all remaining voxels in current row
	  for ( int pX = startX; pX<endX; ++pX, ++r ) 
	    {
	    (pFloating = rowStart) += transformedAxesX[pX];
	    
	    // Continue metric computation.
	    This->ContinueMetric( metricData, r, pFloating );
	    }
	  r += (dimsX-endX);
	  } 
	else
	  {
	  r += dimsX;
	  }
	}
      
      r += (dimsY-endY) * dimsX;
      } 
    else
      {
      r += dimsY * dimsX;
      }
    }

#ifdef CMTK_BUILD_SMP
  This->m_MetricDataMutex.Lock();
#endif
  This->m_MetricData += metricData;
#ifdef CMTK_BUILD_SMP
  This->m_MetricDataMutex.Unlock();
#endif
}

template<class TMultiChannelMetricFunctional>
bool
AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional>
::ClipZ ( const VolumeClipping& clipper, const Vector3D& origin, int& start, int& end ) const
{
  // perform clipping
  Types::Coordinate fromFactor, toFactor;
  if (! clipper.ClipZ( fromFactor, toFactor, origin ) )
    return false;
  
  // there is an intersection: Look up the corresponding grid indices
  start = static_cast<int>( (this->m_ReferenceDims[2]-1)*fromFactor );
  end = 1+std::min( (int)(this->m_ReferenceDims[2]-1), (int)(1 + ((this->m_ReferenceDims[2]-1)*toFactor) ) );
  
  // finally, apply cropping boundaries of the reference volume
  start = std::max<int>( start, this->m_ReferenceCropRegion.From()[2] );
  end = std::min<int>( end, this->m_ReferenceCropRegion.To()[2] );
  
  // return true iff index range is non-empty.
  return (start < end );
}

template<class TMultiChannelMetricFunctional>
bool
AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional>
::ClipX ( const VolumeClipping& clipper, const Vector3D& origin, int& start, int& end ) const
{
  // perform clipping
  Types::Coordinate fromFactor, toFactor;
  if ( ! clipper.ClipX( fromFactor, toFactor, origin, 0, 2, false, true ) )
    return false;
  
  fromFactor = std::min<Types::Coordinate>( 1.0, fromFactor );
  
  // there is an intersection: Look up the corresponding grid indices
  start = std::max( 0, (int)((this->m_ReferenceDims[0]-1)*fromFactor)-1 );
  while ( ( start*this->m_ReferenceChannels[0]->m_Delta[0] < fromFactor*this->m_ReferenceSize[0]) && ( start < this->m_ReferenceDims[0] ) ) 
    ++start;
  
  if ( (toFactor > 1.0) || (start == this->m_ReferenceDims[0]) ) 
    {
    end = this->m_ReferenceDims[0];
    } 
  else
    {
    end = std::min( this->m_ReferenceDims[0]-2, (int)(1 + (this->m_ReferenceDims[0]-1)*toFactor));
    while ( end*this->m_ReferenceChannels[0]->m_Delta[0] > toFactor*this->m_ReferenceSize[0] ) // 'if' not sufficient!	
      --end;
    ++end; // otherwise end=1+min(...) and ...[0][end-1] above!!
    }
  
  // finally, apply cropping boundaries of the reference volume
  start = std::max<int>( start, this->m_ReferenceCropRegion.From()[0] );
  end = std::min<int>( end, this->m_ReferenceCropRegion.To()[0] );
  
  // return true iff index range is non-empty.
  return (start < end );
}

template<class TMultiChannelMetricFunctional>
bool
AffineMultiChannelRegistrationFunctional<TMultiChannelMetricFunctional>
::ClipY ( const VolumeClipping& clipper, const Vector3D& origin, int& start, int& end ) const
{
  // perform clipping
  Types::Coordinate fromFactor, toFactor;
  if ( !clipper.ClipY( fromFactor, toFactor, origin ) )
    return false;
  
  // there is an intersection: Look up the corresponding grid indices
  start = static_cast<int>( (this->m_ReferenceDims[1]-1)*fromFactor );
  
  if ( toFactor > 1.0 ) 
    {
    end = this->m_ReferenceDims[1];
    } 
  else
    {
    end = 1+std::min( this->m_ReferenceDims[1]-1, (int)(1+(this->m_ReferenceDims[1]-1)*toFactor ) );
    }
  // finally, apply cropping boundaries of the reference volume
  start = std::max<int>( start, this->m_ReferenceCropRegion.From()[1] );
  end = std::min<int>( end, this->m_ReferenceCropRegion.To()[1] );
  
  // return true iff index range is non-empty.
  return (start < end );
}

} // namespace cmtk
