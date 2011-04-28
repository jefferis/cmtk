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
#include <System/cmtkConsole.h>

#include <Registration/cmtkReformatVolume.h>

#include <Base/cmtkVector.h>
#include <Base/cmtkAffineXform.h>
#include <Base/cmtkWarpXform.h>
#include <Base/cmtkSplineWarpXform.h>
#include <Base/cmtkSplineWarpXformUniformVolume.h>
#include <Base/cmtkAnatomicalOrientation.h>
#include <Base/cmtkTypedArray.h>
#include <Base/cmtkTemplateArray.h>
#include <Base/cmtkMathUtil.h>
#include <Base/cmtkUniformVolumeInterpolator.h>
#include <Base/cmtkSincInterpolator.h>
#include <Base/cmtkLinearInterpolator.h>
#include <Base/cmtkNearestNeighborInterpolator.h>
#include <Base/cmtkCubicInterpolator.h>
#include <Base/cmtkUniformVolumeInterpolatorPartialVolume.h>

#include <cassert>
#include <algorithm>
#include <vector>

namespace
cmtk
{

/** \addtogroup Registration */
//@{
ReformatVolume::ReformatVolume() 
  : m_UserDataType( TYPE_NONE ),
    ReferenceVolume( NULL ),
    FloatingVolume( NULL ),
    m_AffineXform( NULL ),
    m_WarpXform( NULL )
{
  Interpolation = cmtk::Interpolators::LINEAR;
  this->m_UsePaddingValue = false;
}

void
ReformatVolume::SetReferenceVolume
( const UniformVolume::SmartConstPtr& referenceVolume )
{
  this->ReferenceVolume = referenceVolume;
}

void 
ReformatVolume::SetFloatingVolume
( const UniformVolume::SmartConstPtr& floatingVolume )
{
  FloatingVolume = floatingVolume;
}

void
ReformatVolume::SetAffineXform( const AffineXform::SmartPtr& affineXform )
{
  this->m_AffineXform = affineXform;
}

void
ReformatVolume::SetWarpXform( const WarpXform::SmartPtr& warpXform )
{
  this->m_WarpXform = warpXform;
}

const UniformVolume::SmartPtr
ReformatVolume::MakeTargetVolume() const
{
  return UniformVolume::SmartPtr( ReferenceVolume->CloneGrid() );
}

const UniformVolume::SmartPtr
ReformatVolume::PlainReformat()
{
  UniformVolume::SmartPtr targetVolume = this->MakeTargetVolume();

  if ( targetVolume ) 
    {
    Progress::Begin( 0, targetVolume->GetDims()[AXIS_Z], 1, "Volume reformatting" );
    
    TypedArray::SmartPtr targetData( TypedArray::Create( FloatingVolume->GetData()->GetType(), targetVolume->GetNumberOfPixels() ) );
    if ( this->m_UsePaddingValue )
      {
      targetData->SetPaddingValue( this->m_PaddingValue );
      }
    else
      {
      if ( FloatingVolume->GetData()->GetPaddingFlag() )
	{
	targetData->SetPaddingValue( FloatingVolume->GetData()->GetPaddingValue() );
	}
      }
    
    UniformVolumeInterpolatorBase::SmartPtr interpolator( this->CreateInterpolator( this->FloatingVolume ) );
    Vector3D pFlt;
    
    const DataGrid::IndexType dims = targetVolume->GetDims();
    
    size_t offset = 0;
    for ( int pZ = 0; pZ < dims[2]; ++pZ ) 
      {
      Types::DataItem value = 0;
      
      const SplineWarpXform::SmartConstPtr& splineWarp = SplineWarpXform::SmartConstPtr::DynamicCastFrom( this->m_WarpXform );
      if ( splineWarp ) 
	{
	SplineWarpXformUniformVolume xformVolume( *(this->ReferenceVolume), splineWarp );
	
	for ( int pY = 0; pY<dims[1]; ++pY ) 
	  {
	  for ( int pX = 0; pX<dims[0]; ++pX, ++offset ) 
	    {
	    xformVolume.GetTransformedGrid( pFlt, pX, pY, pZ );
	    
	    if ( interpolator->GetDataAt( pFlt, value ) )
	      targetData->Set( value, offset );	      
	    else
	      targetData->SetPaddingAt( offset );
	    }
	  }
	} 
      else
	{
	for ( int pY = 0; pY<dims[1]; ++pY ) 
	  {
	  for ( int pX = 0; pX<dims[0]; ++pX, ++offset ) 
	    {	    
	    pFlt = ReferenceVolume->GetGridLocation( pX, pY, pZ );
	    this->m_AffineXform->ApplyInPlace( pFlt );
	    
	    if ( interpolator->GetDataAt( pFlt, value ) )
	      targetData->Set( value, offset );
	    else
	      targetData->SetPaddingAt( offset );
	    }
	  }
	}
      
      Progress::SetProgress( pZ );
      }
    
    targetVolume->SetData( targetData );
    }
  
  return targetVolume;
}

TypedArray::SmartPtr
ReformatVolume::PlainReformat
( const int plane, TypedArray::SmartPtr& target, const size_t targetOffset )
{
  const DataGrid::IndexType& Dims = ReferenceVolume->GetDims();
  const int DimsX = Dims[0], DimsY = Dims[1];
  const int DataSize = DimsX * DimsY;

  TypedArray::SmartPtr result = target;
  if ( ! result ) 
    {
    result = TypedArray::Create( FloatingVolume->GetData()->GetType(), DataSize );
     
    if ( this->m_UsePaddingValue )
      result->SetPaddingValue( this->m_PaddingValue );
    }
  
  if ( ! result ) return result;
  
  Vector3D pMod;

  Types::DataItem value = 0;
  int offset = targetOffset;

  UniformVolumeInterpolatorBase::SmartPtr interpolator( this->CreateInterpolator( this->FloatingVolume ) );
  
  const SplineWarpXform::SmartConstPtr splineWarp = SplineWarpXform::SmartConstPtr::DynamicCastFrom( this->m_WarpXform );
  if ( splineWarp ) 
    {
    const SplineWarpXformUniformVolume xformVolume( *(this->ReferenceVolume), splineWarp );
    
    for ( int pY = 0; pY<DimsY; ++pY ) 
       {
      for ( int pX = 0; pX<DimsX; ++pX, ++offset ) 
 	{
	xformVolume.GetTransformedGrid( pMod, pX, pY, plane );
 	
	if ( interpolator->GetDataAt( pMod, value ) )
	  result->Set( value, offset );	      
	else
	  result->SetPaddingAt( offset );
 	}
       }
    } 
  else
    {
    if ( ! this->m_AffineXform ) return result;
    
    for ( int pY = 0; pY<DimsY; ++pY ) 
      {
      for ( int pX = 0; pX<DimsX; ++pX, ++offset ) 
	{	    
	pMod = ReferenceVolume->GetGridLocation( pX, pY, plane );
	this->m_AffineXform->ApplyInPlace( pMod );

	if ( interpolator->GetDataAt( pMod, value ) )
	  result->Set( value, offset );
	else
	  result->SetPaddingAt( offset );
	}
      }
    }
  
  return result;
}

UniformVolume* 
ReformatVolume::GetTransformedReference
( Types::Coordinate *const volumeOffset )
{
  UniformVolume* result = NULL;

  const SplineWarpXform* splineXform = dynamic_cast<const SplineWarpXform*>( this->m_WarpXform.GetConstPtr() );
  if ( ! splineXform ) 
    {
    StdErr << "ERROR: ReformatVolume::GetTransformedReference supports spline warp only.\n";
    return NULL;
    }

  Types::Coordinate bbFrom[3], delta[3];
  result = this->CreateTransformedReference( bbFrom, delta, volumeOffset );

  ScalarDataType dtype = ReferenceVolume->GetData()->GetType();
  TypedArray::SmartPtr dataArray( TypedArray::Create( dtype, result->GetNumberOfPixels() ) );

  if ( this->m_UsePaddingValue )
    dataArray->SetPaddingValue( this->m_PaddingValue );

  result->SetData( dataArray );

  GetTransformedReferenceTP params;
  params.thisObject = this;
  params.ThisThreadIndex = 0;
  params.NumberOfThreads = 1;
  params.dims = result->GetDims();
  params.bbFrom = bbFrom;
  params.delta = delta;
  params.splineXform = splineXform;
  params.dataArray = dataArray;

  UniformVolumeInterpolatorBase::SmartPtr interpolator( this->CreateInterpolator( this->FloatingVolume ) );
  params.referenceInterpolator = interpolator;
  
  DataClass dataClass = ReferenceVolume->GetData()->GetDataClass();
  switch ( dataClass ) 
    {
    default:
    case DATACLASS_GREY: 
    {
    GetTransformedReferenceGrey( &params );
    }
    break;
    case DATACLASS_LABEL: 
    {
    GetTransformedReferenceLabel( &params );
    }
    break;
    }
  
  return result;
}

CMTK_THREAD_RETURN_TYPE
ReformatVolume::GetTransformedReferenceGrey( void *const arg )
{
  GetTransformedReferenceTP* params = static_cast<GetTransformedReferenceTP*>( arg );

  TypedArray::SmartPtr dataArray = params->dataArray;
  const SplineWarpXform* splineXform = params->splineXform;
  const UniformVolumeInterpolatorBase* interpolator = params->referenceInterpolator;
  const Types::Coordinate* delta = params->delta;
  const Types::Coordinate* bbFrom = params->bbFrom;
  const DataGrid::IndexType& dims = params->dims;

  const Types::Coordinate minDelta = MathUtil::Min( 3, delta );
  
  Types::DataItem value;

  UniformVolume::CoordinateVectorType u, v;
  v[2] = bbFrom[2];
  size_t offset = 0;
  for ( int cz = 0; cz < dims[2]; ++cz, v[2] += delta[2] ) 
    {
    if ( ! params->ThisThreadIndex ) Progress::SetProgress( cz );
    v[1] = bbFrom[1];
    for ( int cy = 0; cy < dims[1]; ++cy, v[1] += delta[1] ) 
      {
      v[0] = bbFrom[0];
      for ( int cx = 0; cx < dims[0]; ++cx, v[0] += delta[0], ++offset ) 
	{
	u = v;
	const bool success = splineXform->ApplyInverseInPlace( u, 0.1 * minDelta );
	
	if ( success ) 
	  {
	  if ( interpolator->GetDataAt( u, value ) )
	    dataArray->Set( value, offset );
	  else
	    dataArray->SetPaddingAt( offset );
	  }
	}
      }
    }
  
  return CMTK_THREAD_RETURN_VALUE;
}

CMTK_THREAD_RETURN_TYPE
ReformatVolume::GetTransformedReferenceLabel( void *const arg )
{
  GetTransformedReferenceTP* params = static_cast<GetTransformedReferenceTP*>( arg );

  const ReformatVolume* thisObject = params->thisObject;
  TypedArray::SmartPtr dataArray = params->dataArray;
  const SplineWarpXform* splineXform = params->splineXform;
  const Types::Coordinate* delta = params->delta;
  const Types::Coordinate* bbFrom = params->bbFrom;
  const DataGrid::IndexType& dims = params->dims;

  const std::vector<SplineWarpXform::SmartPtr>* xformList = params->xformList;
  const std::vector<UniformVolume::SmartPtr>* volumeList = params->volumeList;

  const Types::Coordinate minDelta = MathUtil::Min( 3, delta );
  
  UniformVolume::CoordinateVectorType u, v, xyz;

  std::vector<ProbeInfo> probe( params->numberOfImages );
  std::vector<Types::Coordinate> labelCount( params->maxLabel+1 );
    
  xyz[2] = bbFrom[2];
  size_t offset = 0;
  for ( int cz = 0; cz < dims[2]; ++cz, xyz[2] += delta[2] ) 
    {
    if ( ! params->ThisThreadIndex ) Progress::SetProgress( cz );
    xyz[1] = bbFrom[1];
    for ( int cy = 0; cy < dims[1]; ++cy, xyz[1] += delta[1] ) 
      {
      xyz[0] = bbFrom[0];
      for ( int cx = 0; cx < dims[0]; ++cx, xyz[0] += delta[0], ++offset ) 
	{
	v = xyz;
	const bool success = splineXform->ApplyInverseInPlace( v, 0.1 * minDelta );
	u = v;
	
	bool valid = false;
	unsigned int toIdx = 0;
	if ( success ) 
	  {
	  if ( params->IncludeReferenceData ) 
	    {
	    valid = thisObject->ReferenceVolume->ProbeNoXform( probe[toIdx], v );
	    if ( valid ) ++toIdx;
	    }
	  
	  for ( unsigned int img = 0; img < params->numberOfImages-1; ++img ) 
	    {
	    v = u;
	    (*xformList)[img]->ApplyInPlace( v );
	    valid = (*volumeList)[img]->ProbeNoXform( probe[toIdx], v );
	    if ( valid ) ++toIdx;	    
	    }
	  }
	if ( toIdx && success ) 
	  {
	  std::fill( labelCount.begin(), labelCount.end(), 0 );
	  for ( unsigned int idx = 0; idx < toIdx; ++idx ) 
	    {
	    for ( unsigned int corner = 0; corner < 8; ++corner ) 
	      {
	      labelCount[static_cast<int>( probe[idx].Values[corner] )] += probe[idx].GetWeight( corner );
	      }
	    }
	  unsigned int winner = 0;
	  Types::Coordinate winnerWeight = labelCount[0];
	  for ( int label = 1; label < params->maxLabel; ++label ) 
	    {
	    if ( labelCount[label] > winnerWeight ) 
	      {
	      winnerWeight = labelCount[label];
	      winner = label;
	      }
	    }
	  dataArray->Set( static_cast<Types::DataItem>( winner ), offset );
	  } 
	else
	  dataArray->SetPaddingAt( offset );
	}
      }
    }
  
  return CMTK_THREAD_RETURN_VALUE;
}

UniformVolume* 
ReformatVolume::GetTransformedReferenceJacobianAvg
( const std::vector<SplineWarpXform::SmartPtr>* xformList,
  Types::Coordinate *const volumeOffset,
  const bool includeReferenceData )
{
  const SplineWarpXform* splineXform = dynamic_cast<const SplineWarpXform*>( this->m_WarpXform.GetConstPtr() );
  if ( ! splineXform ) 
    {
    StdErr << "ERROR: ReformatVolume::GetTransformedReferenceJacobian supports spline warp only.\n";
    return NULL;
    }
  
  // bounding box for reformatted volume.
  Types::Coordinate bbFrom[3], delta[3];
  UniformVolume* result = this->CreateTransformedReference( bbFrom, delta, volumeOffset );
  
  TypedArray::SmartPtr dataArray( TypedArray::Create( TYPE_FLOAT, result->GetNumberOfPixels() ) );

  if ( this->m_UsePaddingValue )
    dataArray->SetPaddingValue( this->m_PaddingValue );

  result->SetData( dataArray );

  const size_t numberOfThreads = Threads::GetNumberOfThreads();
  std::vector<GetTransformedReferenceTP> params( numberOfThreads );

  for ( size_t thr = 0; thr < numberOfThreads; ++thr ) 
    {
    params[thr].thisObject = this;
    params[thr].ThisThreadIndex = thr;
    params[thr].NumberOfThreads = numberOfThreads;
    params[thr].dims = result->GetDims();
    params[thr].bbFrom = bbFrom;
    params[thr].delta = delta;
    params[thr].splineXform = splineXform;
    params[thr].xformList = xformList;
    params[thr].dataArray = dataArray;
    params[thr].avgMode = MODE_MEAN;
    params[thr].IncludeReferenceData = includeReferenceData;
    }
  
  Threads::RunThreads( GetTransformedReferenceJacobianAvgThread, numberOfThreads, &params[0] );
  
  return result;
}

CMTK_THREAD_RETURN_TYPE
ReformatVolume::GetTransformedReferenceJacobianAvgThread
( void *const arg )
{
  GetTransformedReferenceTP* params = static_cast<GetTransformedReferenceTP*>( arg );

  TypedArray::SmartPtr dataArray = params->dataArray;
  const SplineWarpXform* splineXform = params->splineXform;
  const Types::Coordinate* delta = params->delta;
  const Types::Coordinate* bbFrom = params->bbFrom;
  const DataGrid::IndexType& dims = params->dims;

  const std::vector<SplineWarpXform::SmartPtr>* xformList = params->xformList;

  const Types::Coordinate minDelta = MathUtil::Min( 3, delta );
  
  UniformVolume::CoordinateVectorType u, v, xyz;

  const size_t numberOfXforms = xformList->size();
  std::vector<const SplineWarpXform*> xforms( numberOfXforms );
  for ( unsigned int img = 0; img < numberOfXforms; ++img )
    xforms[img] = (*xformList)[img];

  const int czFrom = params->ThisThreadIndex * dims[2] / params->NumberOfThreads;
  const int czTo = std::min<int>( dims[2], (1+params->ThisThreadIndex) * dims[2] / params->NumberOfThreads );

  Vector<Types::Coordinate> values( params->IncludeReferenceData ? numberOfXforms+1 : numberOfXforms );
  const size_t margin = numberOfXforms / 20; // margin for center 90%

  xyz[2] = bbFrom[2] + czFrom * delta[2];
  size_t offset = czFrom * dims[0] * dims[1];

  for ( int cz = czFrom; cz < czTo; ++cz, xyz[2] += delta[2] ) 
    {
    if ( ! params->ThisThreadIndex ) Progress::SetProgress( cz );
    xyz[1] = bbFrom[1];
    for ( int cy = 0; cy < dims[1]; ++cy, xyz[1] += delta[1] ) 
      {
      xyz[0] = bbFrom[0];
      for ( int cx = 0; cx < dims[0]; ++cx, xyz[0] += delta[0], ++offset ) 
	{
	v = xyz;
	const bool success = splineXform->ApplyInverseInPlace( v, 0.1 * minDelta );
	u = v;
	
	if ( success ) 
	  {
	  const Types::Coordinate refJacobian = splineXform->GetGlobalScaling() / splineXform->GetJacobianDeterminant( u );
	  
	  switch ( params->avgMode ) 
	    {
	    case MODE_MEAN: 
	    {
	    // average
			Types::Coordinate sum = params->IncludeReferenceData ? 1.0 : 0.0;
	    for ( unsigned int img = 0; img < numberOfXforms; ++img )
	      sum += xforms[img]->GetJacobianDeterminant( u ) / xforms[img]->GetGlobalScaling();
	    dataArray->Set( static_cast<Types::DataItem>( refJacobian * sum / numberOfXforms ), offset );
	    break;
	    }
	    case MODE_MEDIAN:
	    case MODE_ROBUST90: 
	    {
	    for ( unsigned int img = 0; img < numberOfXforms; ++img )
	      values[img] = xforms[img]->GetJacobianDeterminant( u ) / xforms[img]->GetGlobalScaling();
	    if ( params->IncludeReferenceData )
	      values[numberOfXforms] = 1.0;
	    
	    values.Sort();
	    
	    if ( params->avgMode == MODE_MEDIAN ) 
	      {
	      // median
	      if ( numberOfXforms & 1 )
		dataArray->Set( static_cast<Types::DataItem>( refJacobian * values[numberOfXforms/2+1] ), offset );
	      else
		dataArray->Set( static_cast<Types::DataItem>( 0.5 * refJacobian * (values[numberOfXforms/2]+values[numberOfXforms/2+1]) ), offset );
	      } 
	    else
	      {
	      // robust average of center 90% average
	      Types::Coordinate sum = 0;
	      for ( unsigned int img = margin; img < numberOfXforms - margin; ++img )
		sum += values[img];
		  dataArray->Set( static_cast<Types::DataItem>( refJacobian * sum / ( numberOfXforms - 2 * margin ) ), offset );
	      }
	    break;
	    }
	    }
	  } 
	else
	  {
	  dataArray->SetPaddingAt( offset );
	  }
	}
      }
    }
  
  return CMTK_THREAD_RETURN_VALUE;
}

UniformVolume*
ReformatVolume::CreateTransformedReference
( Types::Coordinate *const bbFrom, Types::Coordinate *const delta, 
  Types::Coordinate *const volumeOffset )
{
  // default: no offset output desired, 
  // that means, reformat to original image domain
  UniformVolume::CoordinateVectorType bbTo;
  for ( unsigned int axis = 0; axis < 3; ++axis ) 
    {
    bbFrom[axis] = 0;
    bbTo[axis] = this->ReferenceVolume->Size[axis];
    }
  
  if ( volumeOffset ) 
    {
    Vector3D u, v;
    for ( unsigned int z = 0; z < 2; ++z ) 
      {
      if ( z )
	u[AXIS_Z] = this->ReferenceVolume->Size[AXIS_Z];
      else
	u[AXIS_Z] = 0;
      
      for ( unsigned int y = 0; y < 2; ++y ) 
	{
	if ( y )
	  u[AXIS_Y] = this->ReferenceVolume->Size[AXIS_Y];
	else
	  u[AXIS_Y] = 0;
	for ( unsigned int x = 0; x < 2; ++x ) 
	  {
	  if ( x )
	    u[AXIS_X] = this->ReferenceVolume->Size[AXIS_X];
	  else
	    u[AXIS_X] = 0;
	  
	  v = this->m_WarpXform->Apply( u );
	  for ( unsigned int axis = 0; axis < 3; ++axis ) 
	    {
	    bbFrom[axis] = std::min( bbFrom[axis], v[axis] );
	    bbTo[axis] = std::max( bbTo[axis], v[axis] );
	    }
	  }
	}
      }
    
    // report new volume offset, if so desired by the caller
    for ( unsigned int axis = 0; axis < 3; ++axis )
      volumeOffset[axis] = bbFrom[axis];
    }
  
  UniformVolume::IndexType dims;
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    delta[dim] = this->ReferenceVolume->m_Delta[dim];
    bbTo[dim] -= bbFrom[dim];
    dims[dim] = 1 + static_cast<int>( bbTo[dim] / delta[dim] );
    }
  
  return new UniformVolume( dims, bbTo );
}

UniformVolumeInterpolatorBase::SmartPtr
ReformatVolume::CreateInterpolator
( const cmtk::Interpolators::InterpolationEnum interpolation, const UniformVolume::SmartConstPtr& volume )
{
  switch ( interpolation )
    {
    default:
    case cmtk::Interpolators::LINEAR:
    {
    typedef UniformVolumeInterpolator<cmtk::Interpolators::Linear> TInterpolator;
    return TInterpolator::SmartPtr( new TInterpolator( *volume ) );
    }
    case cmtk::Interpolators::NEAREST_NEIGHBOR:
    {
    typedef UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor> TInterpolator;
    return TInterpolator::SmartPtr( new TInterpolator( *volume ) );
    }
    case cmtk::Interpolators::CUBIC:
    {
    typedef UniformVolumeInterpolator<cmtk::Interpolators::Cubic> TInterpolator;
    return TInterpolator::SmartPtr( new TInterpolator( *volume ) );
    }
    case cmtk::Interpolators::COSINE_SINC:
    {
    typedef UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<> > TInterpolator;
    return TInterpolator::SmartPtr( new TInterpolator( *volume ) );
    }
    case cmtk::Interpolators::PARTIALVOLUME:
    {
    typedef UniformVolumeInterpolatorPartialVolume TInterpolator;
    return TInterpolator::SmartPtr( new TInterpolator( *volume ) );
    }
    }  
  return UniformVolumeInterpolatorBase::SmartPtr( NULL );
}

UniformVolumeInterpolatorBase::SmartPtr
ReformatVolume::CreateInterpolator
( const UniformVolume::SmartConstPtr& volume )
{
  return Self::CreateInterpolator( this->Interpolation, volume );
}

} // end namespace cmtk

#ifdef CMTK_BUILD_MPI
#  include "cmtkReformatVolumeMPI.txx"
#else
#  include "cmtkReformatVolume.txx"
#endif
