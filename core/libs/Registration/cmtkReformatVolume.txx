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
//  $Revision: 5806 $
//
//  $LastChangedDate: 2009-05-29 13:36:00 -0700 (Fri, 29 May 2009) $
//
//  $LastChangedBy: torsten $
//
*/

namespace
cmtk
{

/** \addtogroup Registration */
//@{

UniformVolume* 
ReformatVolume::GetTransformedReference
( const std::vector<SplineWarpXform::SmartPtr>* xformList,
  std::vector<UniformVolume::SmartPtr>* volumeList,
  Types::Coordinate *const volumeOffset,
  const bool includeReferenceData )
{
  UniformVolume* result = NULL;
  unsigned int numberOfImages = 0;

  std::vector<const UniformVolumeInterpolatorBase*> interpolatorList;
  interpolatorList.push_back( this->CreateInterpolator( this->ReferenceVolume ) );
  if ( volumeList )
    {
    numberOfImages = 1 + volumeList->size();
    for ( unsigned int img = 0; img < numberOfImages-1; ++img ) 
      {
      interpolatorList.push_back( this->CreateInterpolator( (*volumeList)[img] ) );
      }
    }
  
  const SplineWarpXform* splineXform = dynamic_cast<const SplineWarpXform*>( WarpXform.GetPtr() );
  if ( ! splineXform ) 
    {
    StdErr << "ERROR: ReformatVolume::GetTransformedReference supports spline warp only.\n";
    return NULL;
    }
  
  DataClass dataClass = ReferenceVolume->GetData()->GetDataClass();
  int maxLabel = 0;
  if ( dataClass == DATACLASS_LABEL ) 
    {
    Types::DataItem min, max;
    ReferenceVolume->GetData()->GetRange( min, max );
    maxLabel = static_cast<int>( max );
    
    for ( unsigned int img = 0; img < numberOfImages-1; ++img ) 
      {
      (*volumeList)[img]->GetData()->GetRange( min, max );
      maxLabel = std::max( maxLabel, static_cast<int>( max ) );
      }
    }
  
  // bounding box for reformatted volume.
  Types::Coordinate bbFrom[3], delta[3];
  result = this->CreateTransformedReference( bbFrom, delta, volumeOffset );

  const ScalarDataType dtype = (this->m_UserDataType != TYPE_NONE) ? this->m_UserDataType : ReferenceVolume->GetData()->GetType();
  TypedArray::SmartPtr dataArray( TypedArray::Create( dtype, result->GetNumberOfPixels() ) );
  if ( this->UsePaddingValue )
    dataArray->SetPaddingValue( this->PaddingValue );
  result->SetData( dataArray );

  Progress::SetTotalSteps( result->GetNumberOfPixels() );

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  Array<GetTransformedReferenceTP> params( numberOfThreads );

  for ( size_t thr = 0; thr < numberOfThreads; ++thr ) 
    {
    params[thr].thisObject = this;
    params[thr].ThisThreadIndex = thr;
    params[thr].NumberOfThreads = numberOfThreads;
    params[thr].dims = result->GetDims();
    params[thr].bbFrom = bbFrom;
    params[thr].delta = delta;
    params[thr].splineXform = splineXform;
    params[thr].numberOfImages = numberOfImages;
    params[thr].xformList = xformList;
    params[thr].interpolatorList = &interpolatorList;
    params[thr].dataArray = dataArray;
    params[thr].maxLabel = maxLabel;
    params[thr].IncludeReferenceData = includeReferenceData;
    }
  
  switch ( dataClass ) 
    {
    default:
    case DATACLASS_GREY: 
    {
    if ( xformList && xformList->size() )
      Threads::RunThreads( GetTransformedReferenceGreyAvg, numberOfThreads, &params[0] );
    else
      Threads::RunThreads( GetTransformedReferenceGrey, numberOfThreads, &params[0] );
    }
    break;
    case DATACLASS_LABEL: 
    {
    Threads::RunThreads( GetTransformedReferenceLabel, numberOfThreads, &params[0] );
    }
    break;
    }
  Progress::Done();  
  
  return result;
}

CMTK_THREAD_RETURN_TYPE
ReformatVolume::GetTransformedReferenceGreyAvg( void *const arg )
{
  GetTransformedReferenceTP* params = static_cast<GetTransformedReferenceTP*>( arg );

  TypedArray::SmartPtr dataArray = params->dataArray;
  const SplineWarpXform* splineXform = params->splineXform;
  const Types::Coordinate* delta = params->delta;
  const Types::Coordinate* bbFrom = params->bbFrom;
  const int* dims = params->dims;

  const std::vector<SplineWarpXform::SmartPtr>* xformList = params->xformList;
  const std::vector<const UniformVolumeInterpolatorBase*>* interpolatorList = params->interpolatorList;

  Types::Coordinate minDelta = std::min( delta[0], std::min( delta[1], delta[2] ) );
  
  Array<Types::DataItem> value( params->numberOfImages );
  
  Array<const SplineWarpXform*> xforms( params->numberOfImages-1 );

  for ( unsigned int img = 0; img < params->numberOfImages-1; ++img ) 
    {
    xforms[img] = (*xformList)[img];
    }

  int cx = params->ThisThreadIndex % dims[0];
  int cz = (params->ThisThreadIndex / dims[0]) % dims[1];
  int cy = params->ThisThreadIndex / (dims[0] * dims[1]) ;

  Types::Coordinate x = bbFrom[0] + cx * delta[0];
  Types::Coordinate y = bbFrom[1] + cy * delta[1];
  Types::Coordinate z = bbFrom[2] + cz * delta[2];
  
  const size_t numberOfPixels = dims[0] * dims[1] * dims[2];
  const size_t statusUpdateIncrement = numberOfPixels / 100;

  Vector3D u, v;
  bool success = false;
  for ( size_t offset = params->ThisThreadIndex; offset < numberOfPixels; offset += params->NumberOfThreads ) 
    {
    if ( ! params->ThisThreadIndex && ! (offset % statusUpdateIncrement) ) 
      Progress::SetProgress( offset );

    v.Set( x, y, z );
    if ( cx && success ) 
      {
      // use previous (successful) inverse as starting point inside row
      success = splineXform->ApplyInverseInPlaceWithInitial( v, u, 0.1 * minDelta );
      } 
    else
      {
      success = splineXform->ApplyInverseInPlace( v, 0.1 * minDelta );
      }
    u = v;
    
    unsigned int toIdx = 0;
    if ( success ) 
      {
      if ( params->IncludeReferenceData ) 
	{
	if ( (*interpolatorList)[0]->GetDataAt( v, value[toIdx] ) )
	  ++toIdx;
	}
      
      for ( unsigned int img = 0; img < params->numberOfImages-1; ++img ) 
	{
	v = u;
	xforms[img]->ApplyInPlace( v );
	
	if ( (*interpolatorList)[img+1]->GetDataAt( v, value[toIdx] ) )
	  ++toIdx;
	}
      }
    if ( toIdx && success ) 
      {
      Types::DataItem avg = value[0];
      for ( unsigned int idx = 1; idx < toIdx; ++idx )
	avg += value[idx];
      dataArray->Set( avg / toIdx, offset );
      } 
    else
      dataArray->SetPaddingAt( offset );

    cx += params->NumberOfThreads;
    if ( cx >= dims[0] )
      {
      cy += cx / dims[0];
      cx %= dims[0];

      if ( cy >= dims[1] )
	{
	cz += cy / dims[1];
	cy %= dims[1];
	z = bbFrom[2] + cz * delta[2];
	}
      y = bbFrom[1] + cy * delta[1];
      }
    x = bbFrom[0] + cx * delta[0];
    }

  return CMTK_THREAD_RETURN_VALUE;
}

} // namespace cmtk
