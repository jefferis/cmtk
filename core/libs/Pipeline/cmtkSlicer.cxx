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

#include <cmtkSlicer.h>

#include <cmtkAffineXform.h>
#include <cmtkUniformVolumeInterpolator.h>
#include <cmtkLinearInterpolator.h>
#include <cmtkCubicInterpolator.h>
#include <cmtkSincInterpolator.h>
#include <cmtkNearestNeighborInterpolator.h>
#include <cmtkUniformVolumeInterpolatorPartialVolume.h>

namespace 
cmtk
{

/** \addtogroup Pipeline */
//@{

Slicer::Slicer()
{
  ApplyWarp = false;
  this->m_Plane = NULL;
  InterpolationMode = cmtk::Interpolators::LINEAR;
  ZoomFactor = 1;
  TempImage = NULL;
}

Slicer::~Slicer ()
{
  if ( this->m_Plane ) 
    this->m_Plane->Delete();
  if ( TempImage )
    TempImage->Delete();
}

void
Slicer::SetPlane( Plane *const plane ) 
{
  this->ReplaceObject( this->m_Plane, plane );
}

long
Slicer::Update()
{
  this->CheckInputForUpdate( this->m_Plane );
  this->CheckInputForUpdate( Input );
  return this->ExecuteIfNecessary();
}

void
Slicer::Execute()
{
  if ( Input == NULL ) return;
  if ( this->m_Plane == NULL ) return;
  UniformVolume::SmartPtr volume = this->Input->GetVolume();
  if ( ! volume ) return;

  WarpXform::SmartPtr warpXform = Input->GetWarpXform();
  
  Image *image = NULL;

  if ( ZoomFactor == 1 )
    image = this->GetOutput();
  else {
    if ( ! TempImage ) {
      TempImage = Image::New();
    }
    image = TempImage;
  }
  image->CopyStructure( this->m_Plane );

  image->SetDataType( volume->GetData()->GetType() );
  TypedArray::SmartPtr data = image->GetData();

  Vector3D p( this->m_Plane->GetOrigin() );
  
  Vector3D dirX = Vector3D( this->m_Plane->GetDirectionX() );
  Vector3D dirY = Vector3D( this->m_Plane->GetDirectionY() );

  dirX *= (this->m_Plane->GetSpacing()[0] / dirX.RootSumOfSquares());
  dirY *= (this->m_Plane->GetSpacing()[1] / dirY.RootSumOfSquares());

  const unsigned int *dims = this->m_Plane->GetDims();

  cmtk::UniformVolumeInterpolatorBase::SmartPtr interpolator;
  switch ( this->InterpolationMode )
    {
    default:
    case cmtk::Interpolators::LINEAR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::NEAREST_NEIGHBOR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor> TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::CUBIC:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Cubic> TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::COSINE_SINC:
    {
    typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<> > TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::PARTIALVOLUME:
    {
    typedef cmtk::UniformVolumeInterpolatorPartialVolume TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    }

  if ( (!ApplyWarp) || ( !warpXform) ) 
    {
    AffineXform::SmartPtr affineXform( NULL );
    AffineXform::SmartPtr originalXform = Input->GetAffineXform();
    if ( originalXform ) affineXform = originalXform;
    
    if ( affineXform ) 
      {
      affineXform->ApplyInPlace( dirX += p );
      affineXform->ApplyInPlace( dirY += p );
      affineXform->ApplyInPlace( p );
      dirX -= p;
      dirY -= p;
      }
    
    size_t index = 0;
    Types::DataItem value;
    Vector3D rowOffset;
    for ( unsigned int y = 0; y<dims[1]; ++y, p += dirY ) 
      {
      rowOffset = p;
      for ( unsigned int x = 0; x<dims[0]; ++x, p += dirX, ++index ) 
	{
	if ( interpolator->GetDataAt( p, value ) )
	  data->Set( value, index );
	else
	  data->SetPaddingAt( index );
	}
      p = rowOffset;
      }
    }
  else 
    {
    SplineWarpXform::SmartPtr splineWarpXform = SplineWarpXform::SmartPtr::DynamicCastFrom( warpXform );
    if ( splineWarpXform ) 
      {
      this->ExecuteSplineWarp( data, splineWarpXform, dims, p, dirX, dirY );
      }
    else
      {
      fputs( "Coding error: Unsupported WarpXform descendant encountered or dynamic_cast() failed.\n", stderr );
      }
    }

  /// If we are zooming, generate final output by 2-D interpolation.
  if ( ZoomFactor != 1 ) 
    {
    Image *output = this->GetOutput();
    
    unsigned int dims[2];
    image->GetDims( dims );
    
    dims[0] = static_cast<unsigned int>( dims[0] * ZoomFactor );
    dims[1] = static_cast<unsigned int>( dims[1] * ZoomFactor );
    
    Types::Coordinate spacing[2];
    image->GetSpacing( spacing );
    spacing[0] /= ZoomFactor;
    spacing[1] /= ZoomFactor;

    output->SetDims( dims );
    output->SetSpacing( spacing );
    output->SetDataType( image->GetDataType() );

    TypedArray::SmartPtr outData = output->GetData();

    unsigned int offset = 0;
    for ( unsigned int y = 0; y < dims[1]; ++y )
      for ( unsigned int x = 0; x < dims[0]; ++x, ++offset ) 
	{
	outData->Set( image->GetDataAt( x * spacing[0], y * spacing[1] ), offset );
	}
    }
  
  this->UpdateExecuteTime();
}

void
Slicer::ExecuteSplineWarp
( TypedArray::SmartPtr& data, const SplineWarpXform* warpXform,
  const unsigned int* dims, const Vector3D& offset,
  const Vector3D& dX, const Vector3D& dY )
{
  UniformVolume::SmartPtr volume = this->Input->GetVolume();
  if ( ! volume) return;

#ifdef DEBUG
  puts( "Slicer: Applying spline warp." );
#endif
  unsigned int index = 0;
  Types::DataItem value;
  Vector3D pDeformed, rowOffset;
  Vector3D p( offset );

  cmtk::UniformVolumeInterpolatorBase::SmartPtr interpolator;
  switch ( this->InterpolationMode )
    {
    default:
    case cmtk::Interpolators::LINEAR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Linear> TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::NEAREST_NEIGHBOR:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::NearestNeighbor> TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::CUBIC:
    {
    typedef cmtk::UniformVolumeInterpolator<cmtk::Interpolators::Cubic> TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::COSINE_SINC:
    {
    typedef cmtk::UniformVolumeInterpolator< cmtk::Interpolators::CosineSinc<> > TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    case cmtk::Interpolators::PARTIALVOLUME:
    {
    typedef cmtk::UniformVolumeInterpolatorPartialVolume TInterpolator;
    interpolator = cmtk::UniformVolumeInterpolatorBase::SmartPtr( new TInterpolator( volume ) );
    break;
    }
    }

  for ( unsigned int y = 0; y<dims[1]; ++y, p += dY ) 
    {
    rowOffset = p;
    for ( unsigned int x = 0; x<dims[0]; ++x, p += dX, ++index ) 
      {
      pDeformed = p;
      warpXform->ApplyInPlace( pDeformed );
      if ( interpolator->GetDataAt( pDeformed, value ) )
	data->Set( value, index );
      else
	data->SetPaddingAt( index );
      }
    p = rowOffset;
    }
}

} // namespace cmtk
