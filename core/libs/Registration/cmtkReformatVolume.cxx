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

#include <cmtkReformatVolume.h>

#include <cmtkVector.h>
#include <cmtkAffineXform.h>
#include <cmtkWarpXform.h>
#include <cmtkSplineWarpXform.h>
#include <cmtkAnatomicalOrientation.h>
#include <cmtkTypedArray.h>
#include <cmtkTemplateArray.h>
#include <cmtkMathUtil.h>

#include <cmtkProgress.h>
#include <cmtkConsole.h>

#include <assert.h>
#include <algorithm>
#include <vector>

#include <cmtkUniformVolumeInterpolator.h>
#include <cmtkSincInterpolator.h>
#include <cmtkLinearInterpolator.h>
#include <cmtkCubicInterpolator.h>
#include <cmtkUniformVolumeInterpolatorPartialVolume.h>

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
    m_WarpXform( NULL ),
    CheckerboardMode( false )
{
  Interpolation = cmtk::Interpolators::LINEAR;

  RescaleReference = 0;
  RescaleOffset = 0;
  RescaleSlope = 1;

  this->m_UsePaddingValue = false;
}

void
ReformatVolume::SetReferenceVolume
( const UniformVolume::SmartPtr& referenceVolume )
{
  ReferenceVolume = referenceVolume;
  if ( ReferenceVolume && ReferenceVolume->GetData() ) 
    {
    Types::DataItem min, max;
    ReferenceVolume->GetData()->GetRange( min, max );
    MaximumValue = std::max( max, MaximumValue );
    }			       
}

void 
ReformatVolume::SetFloatingVolume
( const UniformVolume::SmartPtr& floatingVolume )
{
  FloatingVolume = floatingVolume;
  if ( FloatingVolume && FloatingVolume->GetData() ) 
    {
    Types::DataItem min, max;
    FloatingVolume->GetData()->GetRange( min, max );
    MaximumValue = std::max( max, MaximumValue );
    }			       
}

void
ReformatVolume::SetAffineXform( const AffineXform::SmartPtr& affineXform )
{
  this->m_AffineXform = affineXform;
}

void
ReformatVolume::SetWarpXform( WarpXform::SmartPtr& warpXform )
{
  this->m_WarpXform = warpXform;
}

void
ReformatVolume::SetRescale
( const int rescaleReference, const Types::DataItem rescaleOffset, const Types::DataItem rescaleSlope )
{
  RescaleReference = rescaleReference;
  RescaleOffset = rescaleOffset;
  RescaleSlope = rescaleSlope;
  this->SetRescale();
}

UniformVolume* 
ReformatVolume::MakeTargetVolume() const
{
  return ReferenceVolume->CloneGrid();
}

UniformVolume* 
ReformatVolume::PlainReformat()
{
  UniformVolume* targetVolume = this->MakeTargetVolume();

  if ( targetVolume ) 
    {
    size_t planeSize = targetVolume->GetDims(AXIS_X) * targetVolume->GetDims(AXIS_Y);
    
    Progress::SetTotalSteps( targetVolume->GetDims( AXIS_Z ) );
    
    TypedArray::SmartPtr targetData( TypedArray::Create( FloatingVolume->GetData()->GetType(), targetVolume->GetNumberOfPixels() ) );
    if ( this->m_UsePaddingValue )
      targetData->SetPaddingValue( this->m_PaddingValue );
    
    size_t planeOffset = 0;
    for ( int plane = 0; plane < targetVolume->GetDims(AXIS_Z); ++plane, planeOffset += planeSize ) 
      {
      this->PlainReformat( plane, targetData.GetPtr(), planeOffset );
      Progress::SetProgress( plane );
      }
    
    targetVolume->SetData( targetData );
    }
  
  return targetVolume;
}

TypedArray* 
ReformatVolume::PlainReformat
( const int plane, TypedArray *const target, const size_t targetOffset )
{
  const int *Dims = ReferenceVolume->GetDims();
  const int DimsX = Dims[0], DimsY = Dims[1];
  const int DataSize = DimsX * DimsY;

  TypedArray *result = target;
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
  
  if ( this->m_WarpXform ) 
    {
    this->m_WarpXform->RegisterVolume( ReferenceVolume );
    
    const SplineWarpXform::SmartPtr& splineWarp = SplineWarpXform::SmartPtr::DynamicCastFrom( this->m_WarpXform );
    
    for ( int pY = 0; pY<DimsY; ++pY ) 
      {
      for ( int pX = 0; pX<DimsX; ++pX, ++offset ) 
	{
	splineWarp->GetTransformedGrid( pMod, pX, pY, plane );
	
	if ( interpolator->GetDataAt( pMod, value ) )
	  result->Set( value, offset );	      
	else
	  if ( CheckerboardMode )
	    result->Set( ((pX>>3)%2)==((pY>>3)%2) ? 0 : MaximumValue, offset );
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
	ReferenceVolume->GetGridLocation( pMod, pX, pY, plane );
	this->m_AffineXform->ApplyInPlace( pMod );

	if ( interpolator->GetDataAt( pMod, value ) )
	  result->Set( value, offset );
	else
	  if ( CheckerboardMode )
	    result->Set( ((pX>>3)%2)==((pY>>3)%2) ? 0 : MaximumValue, offset );
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

  const SplineWarpXform* splineXform = dynamic_cast<const SplineWarpXform*>( this->m_WarpXform.GetPtr() );
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

  Progress::SetTotalSteps( result->GetNumberOfPixels() );
    
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
  Progress::Done();  
  
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
  const int* dims = params->dims;

  const Types::Coordinate minDelta = MathUtil::Min( 3, delta );
  
  Vector3D u, v;
  Types::Coordinate x, y, z;
  bool success = false;

  Types::DataItem value;
  z = bbFrom[2];
  size_t offset = 0;
  for ( int cz = 0; cz < dims[2]; ++cz, z += delta[2] ) 
    {
    if ( ! params->ThisThreadIndex ) Progress::SetProgress( cz );
    y = bbFrom[1];
    for ( int cy = 0; cy < dims[1]; ++cy, y += delta[1] ) 
      {
      x = bbFrom[0];
      for ( int cx = 0; cx < dims[0]; ++cx, x += delta[0], ++offset ) 
	{
	v.Set( x, y, z );
	if ( cx && success ) 
	  // use previous (successful) inverse as starting point inside row
	  success = splineXform->ApplyInverseInPlaceWithInitial( v, u, 0.1 * minDelta );
	else 
	  {
	  success = splineXform->ApplyInverseInPlace( v, 0.1 * minDelta );
	  }
	u = v;
	
	if ( success ) 
	  {
	  if ( interpolator->GetDataAt( v, value ) )
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
  GetTransformedReferenceLabelTP* params = static_cast<GetTransformedReferenceLabelTP*>( arg );

  const ReformatVolume* thisObject = params->thisObject;
  TypedArray::SmartPtr dataArray = params->dataArray;
  const SplineWarpXform* splineXform = params->splineXform;
  const Types::Coordinate* delta = params->delta;
  const Types::Coordinate* bbFrom = params->bbFrom;
  const int* dims = params->dims;

  const std::vector<SplineWarpXform::SmartPtr>* xformList = params->xformList;
  const std::vector<UniformVolume::SmartPtr>* volumeList = params->volumeList;

  const Types::Coordinate minDelta = MathUtil::Min( 3, delta );
  
  Vector3D u, v;
  Types::Coordinate x, y, z;
  bool success = false;

  std::vector<ProbeInfo> probe( params->numberOfImages );
  std::vector<Types::Coordinate> labelCount( params->maxLabel+1 );
    
  z = bbFrom[2];
  size_t offset = 0;
  for ( int cz = 0; cz < dims[2]; ++cz, z += delta[2] ) 
    {
    if ( ! params->ThisThreadIndex ) Progress::SetProgress( cz );
    y = bbFrom[1];
    for ( int cy = 0; cy < dims[1]; ++cy, y += delta[1] ) 
      {
      x = bbFrom[0];
      for ( int cx = 0; cx < dims[0]; ++cx, x += delta[0], ++offset ) 
	{
	v.Set( x, y, z );
	if ( cx && success ) 
	  // use previous (successful) inverse as starting point inside row
	  success = splineXform->ApplyInverseInPlaceWithInitial( v, u, 0.1 * minDelta );
	else
	  {
	  success = splineXform->ApplyInverseInPlace( v, 0.1 * minDelta );
	  }
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
  const SplineWarpXform* splineXform = dynamic_cast<const SplineWarpXform*>( this->m_WarpXform.GetPtr() );
  if ( ! splineXform ) 
    {
    StdErr
      << "ERROR: ReformatVolume::GetTransformedReferenceJacobian supports"
      << "spline warp only.\n";
    return NULL;
    }

  // bounding box for reformatted volume.
  Types::Coordinate bbFrom[3], delta[3];
  UniformVolume* result = this->CreateTransformedReference( bbFrom, delta, volumeOffset );
  
  TypedArray::SmartPtr dataArray( TypedArray::Create( TYPE_FLOAT, result->GetNumberOfPixels() ) );

  if ( this->m_UsePaddingValue )
    dataArray->SetPaddingValue( this->m_PaddingValue );

  result->SetData( dataArray );

  Progress::SetTotalSteps( result->GetDims( AXIS_Z ) );

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
  Progress::Done();  
  
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
  const int* dims = params->dims;

  const std::vector<SplineWarpXform::SmartPtr>* xformList = params->xformList;

  const Types::Coordinate minDelta = MathUtil::Min( 3, delta );
  
  Vector3D u, v;
  Types::Coordinate x, y, z;
  bool success = false;

  const size_t numberOfXforms = xformList->size();
  std::vector<const SplineWarpXform*> xforms( numberOfXforms );
  for ( unsigned int img = 0; img < numberOfXforms; ++img )
    xforms[img] = (*xformList)[img];

  const int czFrom = params->ThisThreadIndex * dims[2] / params->NumberOfThreads;
  const int czTo = std::min<int>( dims[2], (1+params->ThisThreadIndex) * dims[2] / params->NumberOfThreads );

#ifdef __IGNORE__
  StdErr.printf( " Thread #%d from %d to %d\n", params->ThisThreadIndex, czFrom, czTo );
#endif

  Vector<Types::Coordinate> values( params->IncludeReferenceData ? numberOfXforms+1 : numberOfXforms );
  const size_t margin = numberOfXforms / 20; // margin for center 90%

  z = bbFrom[2] + czFrom * delta[2];
  size_t offset = czFrom * dims[0] * dims[1];

  for ( int cz = czFrom; cz < czTo; ++cz, z += delta[2] ) 
    {
    if ( ! params->ThisThreadIndex ) Progress::SetProgress( cz );
    y = bbFrom[1];
    for ( int cy = 0; cy < dims[1]; ++cy, y += delta[1] ) 
      {
      x = bbFrom[0];
      for ( int cx = 0; cx < dims[0]; ++cx, x += delta[0], ++offset ) 
	{
	v.Set( x, y, z );
	if ( cx && success ) 
	  // use previous (successful) inverse as starting point inside row
	  success = splineXform->ApplyInverseInPlaceWithInitial( v, u, 0.1 * minDelta );
	else
	  {
	  success = splineXform->ApplyInverseInPlace( v, 0.1 * minDelta );
	  }
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
	      // robust average of center 90%
	      // average
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
  Types::Coordinate bbTo[3];
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
	u.XYZ[AXIS_Z] = this->ReferenceVolume->Size[AXIS_Z];
      else
	u.XYZ[AXIS_Z] = 0;
      
      for ( unsigned int y = 0; y < 2; ++y ) 
	{
	if ( y )
	  u.XYZ[AXIS_Y] = this->ReferenceVolume->Size[AXIS_Y];
	else
	  u.XYZ[AXIS_Y] = 0;
	for ( unsigned int x = 0; x < 2; ++x ) 
	  {
	  if ( x )
	    u.XYZ[AXIS_X] = this->ReferenceVolume->Size[AXIS_X];
	  else
	    u.XYZ[AXIS_X] = 0;
	  
	  v = this->m_WarpXform->Apply( u );
	  for ( unsigned int axis = 0; axis < 3; ++axis ) 
	    {
	    bbFrom[axis] = std::min( bbFrom[axis], v.XYZ[axis] );
	    bbTo[axis] = std::max( bbTo[axis], v.XYZ[axis] );
	    }
	  }
	}
      }
    
    // report new volume offset, if so desired by the caller
    for ( unsigned int axis = 0; axis < 3; ++axis )
      volumeOffset[axis] = bbFrom[axis];
    }
  
  int dims[3];
  for ( int dim = 0; dim < 3; ++dim ) 
    {
    delta[dim] = this->ReferenceVolume->m_Delta[dim];
    bbTo[dim] -= bbFrom[dim];
    dims[dim] = 1 + static_cast<int>( bbTo[dim] / delta[dim] );
    }
  
  return new UniformVolume( dims, bbTo );
}

UniformVolumeInterpolatorBase*
ReformatVolume::CreateInterpolator
( const cmtk::Interpolators::InterpolationEnum interpolation, const UniformVolume::SmartPtr& volume )
{
  switch ( interpolation )
    {
    default:
    case cmtk::Interpolators::LINEAR:
    {
    typedef UniformVolumeInterpolator<cmtk::Interpolators::Linear> TInterpolator;
    return new TInterpolator( volume );
    break;
    }
    case cmtk::Interpolators::CUBIC:
    {
    typedef UniformVolumeInterpolator<cmtk::Interpolators::Cubic> TInterpolator;
    return new TInterpolator( volume );
    break;
    }
    case cmtk::Interpolators::COSINE_SINC:
    {
    typedef UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<> > TInterpolator;
    return new TInterpolator( volume );
    break;
    }
    case cmtk::Interpolators::PARTIALVOLUME:
    {
    typedef UniformVolumeInterpolatorPartialVolume TInterpolator;
    return new TInterpolator( volume );
    break;
    }
    }  
  return NULL;
}

UniformVolumeInterpolatorBase*
ReformatVolume::CreateInterpolator
( const UniformVolume::SmartPtr& volume )
{
  return Self::CreateInterpolator( this->Interpolation, volume );
}

} // end namespace cmtk

#ifdef CMTK_BUILD_MPI
#  include "cmtkReformatVolumeMPI.txx"
#else
#  include "cmtkReformatVolume.txx"
#endif
