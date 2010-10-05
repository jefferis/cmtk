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

#include "cmtkWarpXform.h"

#include <Base/cmtkMathUtil.h>
#include <Base/cmtkUniformVolume.h>

#include <string.h>
#include <algorithm>

namespace
cmtk
{

/** \addtogroup Base */
//@{

void
WarpXform::InitGrid
( const FixedVector<3,Types::Coordinate>& domain, const DataGrid::IndexType& dims )
{
  this->Domain = domain;
  this->m_Dims = dims;
  std::fill( this->m_Offset.begin(), this->m_Offset.end(), 0 );
  
  NumberOfControlPoints = this->m_Dims[0] * this->m_Dims[1] * this->m_Dims[2];
  this->AllocateParameterVector( 3 * NumberOfControlPoints );
  this->Update();
}

void 
WarpXform::Update( const bool )
{
  nextI = 3;
  nextJ = nextI * this->m_Dims[0];
  nextK = nextJ * this->m_Dims[1];
  nextIJ = nextJ + nextI;
  nextIK = nextK + nextI;
  nextJK = nextK + nextJ;
  nextIJK = nextJK + nextI;
}

void
WarpXform::GetDerivativeLandmarksMSD
( double& lowerMSD, double& upperMSD, const MatchedLandmarkList* ll,
  const unsigned int idx, const Types::Coordinate step )
{
  upperMSD = lowerMSD = 0;

  Types::Coordinate pOld = this->m_Parameters[idx];

  this->m_Parameters[idx] += step;
  MatchedLandmarkList::const_iterator it = ll->begin();
  while ( it != ll->end() ) 
    {
    Self::SpaceVectorType source( (*it)->GetLocation() );
    Self::SpaceVectorType target( (*it)->GetTargetLocation() );
    this->ApplyInPlace( source );
    upperMSD += (source - target).SumOfSquares();
    ++it;
    }
  
  this->m_Parameters[idx] = pOld - step;
  it = ll->begin();
  while ( it != ll->end() ) 
    {
    Self::SpaceVectorType source( (*it)->GetLocation() );
    Self::SpaceVectorType target( (*it)->GetTargetLocation() );
    this->ApplyInPlace( source );
    lowerMSD += (source - target).SumOfSquares();
    ++it;
    }
  this->m_Parameters[idx] = pOld;
  
  upperMSD /= ll->size();
  lowerMSD /= ll->size();
}

Types::Coordinate
WarpXform::GetInverseConsistencyError
( const WarpXform* inverse, const UniformVolume* volume, const DataGrid::RegionType* voi ) const 
{
  Self::SpaceVectorType v, vv;
  Types::Coordinate result = 0.0;
  int count = 0;

  DataGrid::RegionType myVoi;
  const DataGrid::RegionType *pVoi = &myVoi;
  if ( voi ) 
    {
    pVoi = voi;
    } 
  else
    {
    myVoi = volume->GetWholeImageRegion();
    }

  for ( int z = pVoi->From()[2]; z < pVoi->To()[2]; ++z )
    for ( int y = pVoi->From()[1]; y < pVoi->To()[1]; ++y )
      for ( int x = pVoi->From()[0]; x < pVoi->To()[0]; ++x ) 
	{
	v = volume->GetGridLocation( x, y, z );
	vv = v;
	this->ApplyInPlace( vv );
	if ( inverse->InDomain( vv ) ) 
	  {
	  inverse->ApplyInPlace( vv );
	  v -= vv;
	  result += v.RootSumOfSquares();
	  ++count;
	  }
	}
  
  return count ? result / count : 0.0;
}

void
WarpXform::GetDerivativeInverseConsistencyError
( double& lower, double& upper, const WarpXform* inverse,
  const UniformVolume* volume, const DataGrid::RegionType* voi, 
  const unsigned int idx, const Types::Coordinate step )
{
  const Types::Coordinate pOld = this->m_Parameters[idx];
  upper = lower = (-this->GetInverseConsistencyError( inverse, volume, voi ));

  this->m_Parameters[idx] += step;
  upper += this->GetInverseConsistencyError( inverse, volume, voi );

  this->m_Parameters[idx] = pOld - step;
  lower+= this->GetInverseConsistencyError( inverse, volume, voi );

  this->m_Parameters[idx] = pOld;
}

Types::Coordinate 
WarpXform::GetParamStep
( const size_t idx, const Self::SpaceVectorType&, const Types::Coordinate mmStep ) const
{
  if ( this->m_ActiveFlags && ! (*this->m_ActiveFlags)[idx] ) return 0;

  int controlPointIdx = idx / 3;
  unsigned short x =  ( controlPointIdx %  this->m_Dims[0] );
  unsigned short y = ( (controlPointIdx /  this->m_Dims[0]) % this->m_Dims[1] );
  unsigned short z = ( (controlPointIdx /  this->m_Dims[0]) / this->m_Dims[1] );
  
  if ( (x>=this->m_IgnoreEdge) && (x<(this->m_Dims[0]-this->m_IgnoreEdge)) && 
       (y>=this->m_IgnoreEdge) && (y<(this->m_Dims[1]-this->m_IgnoreEdge)) && 
       (z>=this->m_IgnoreEdge) && (z<(this->m_Dims[2]-this->m_IgnoreEdge)) ) 
    {
    return mmStep;
    } 
  else
    {
    return 0;
    }
}

void 
WarpXform::SetParameterActive()
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  this->m_ActiveFlags->Set();
}

void
WarpXform::SetParameterActive
( const size_t index, const bool active )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  this->m_ActiveFlags->Set( index, active );
}

void 
WarpXform::SetParametersActive( const DataGrid::RegionType& )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
}

void
WarpXform::SetParametersActive
( const int axis, const bool active )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  for ( unsigned int idx = (unsigned int)axis; idx < this->m_NumberOfParameters; idx += 3 )
    this->m_ActiveFlags->Set( idx, active );
}

void
WarpXform::SetParametersActive( const char* axes )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  if ( axes ) 
    {
    if ( strchr( axes, 'x' ) || strchr( axes, 'X' ) )
      this->SetParametersActive( AXIS_X );
    if ( strchr( axes, 'y' ) || strchr( axes, 'Y' ) )
      this->SetParametersActive( AXIS_Y );
    if ( strchr( axes, 'z' ) || strchr( axes, 'Z' ) )
      this->SetParametersActive( AXIS_Z );
    }
}

void
WarpXform::SetParameterInactive( const size_t index )
{
  if ( !this->m_ActiveFlags ) 
    {
    this->m_ActiveFlags = BitVector::SmartPtr( new BitVector( this->m_NumberOfParameters, true ) );
    }
  this->m_ActiveFlags->Reset( index );
}

int
WarpXform::GetParameterActive( const size_t index ) const
{
  if ( this->m_ActiveFlags )
    return (*this->m_ActiveFlags)[index];
  else
    return 1;
}

void
WarpXform::DeleteParameterActiveFlags()
{
  this->m_ActiveFlags = BitVector::SmartPtr::Null;
}

void 
WarpXform::SetIncompressibilityMap
( DataGrid::SmartPtr& incompressibility )
{
  IncompressibilityMap = incompressibility;
}

void
WarpXform::ReplaceInitialAffine( const AffineXform* newAffineXform )
{
  AffineXform change;

  // first, get inverse of current initial affine transformation
  if ( this->m_InitialAffineXform ) 
    {
    change = *(this->m_InitialAffineXform->GetInverse());
    }

  // second, append current initial affine transformation
  if ( newAffineXform )
    change.Concat( *newAffineXform );

  // apply effective change to all control points.
  Types::Coordinate *coeff = this->m_Parameters;
  for ( unsigned int idx = 0; idx < NumberOfControlPoints; ++idx, coeff+=3 ) 
    {
    Self::SpaceVectorType p( coeff );
    change.ApplyInPlace( p );
    coeff[0] = p[0];
    coeff[1] = p[1];
    coeff[2] = p[2];
    }

  // Finally, copy new transformation. We want to create a new object here
  // if the current transformation is linked somewhere else.
  if ( newAffineXform )
    {
    this->m_InitialAffineXform = AffineXform::SmartPtr::DynamicCastFrom( newAffineXform->Clone() );
    }
  else
    {
    this->m_InitialAffineXform = AffineXform::SmartPtr( new AffineXform );
    }
  this->m_InitialAffineXform->m_MetaInformation[META_XFORM_FIXED_IMAGE_PATH] = this->m_MetaInformation[META_XFORM_FIXED_IMAGE_PATH];
  this->m_InitialAffineXform->m_MetaInformation[META_XFORM_MOVING_IMAGE_PATH] = this->m_MetaInformation[META_XFORM_MOVING_IMAGE_PATH];
}

void
WarpXform::ConcatAffine( const AffineXform* affineXform )
{
  // apply effective change to all control points.
  Types::Coordinate *coeff = this->m_Parameters;
  for ( unsigned int idx = 0; idx < NumberOfControlPoints; ++idx, coeff+=3 ) 
    {
    Self::SpaceVectorType p( coeff );
    affineXform->ApplyInPlace( p );
    coeff[0] = p[0];
    coeff[1] = p[1];
    coeff[2] = p[2];
    }

  // Finally, generate combined affine transformation. We want to create a new
  // object here if the current transformation is linked somewhere else.
  if ( this->m_InitialAffineXform.GetReferenceCount() != 1 )
    this->m_InitialAffineXform = this->m_InitialAffineXform->Clone();
  this->m_InitialAffineXform->Concat( *affineXform );
}

} // namespace cmtk
